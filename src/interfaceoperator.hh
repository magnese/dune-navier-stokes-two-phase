#ifndef DUNE_FEM_INTERFACEOPERATOR_HH
#define DUNE_FEM_INTERFACEOPERATOR_HH

#define USE_LAGRANGE_QUADRATURE_POINTS 1

#include <dune/geometry/referenceelements.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include "normal.hh"

#include <fstream>
#include <vector>
#include <string>

namespace Dune
{
namespace Fem
{

template<typename LinearOperatorImp>
class InterfaceOperator:public Operator<typename LinearOperatorImp::DomainFunctionType,typename LinearOperatorImp::RangeFunctionType>
{
  public:
  typedef LinearOperatorImp LinearOperatorType;
  typedef typename LinearOperatorType::DomainFunctionType DomainFunctionType;
  typedef typename LinearOperatorType::RangeFunctionType RangeFunctionType;
  typedef DomainFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;
  typedef InterfaceOperator<LinearOperatorType> ThisType;

  explicit InterfaceOperator(const DiscreteSpaceType& space):
    space_(space),op_("interface operator",space_,space_)
  {}

  InterfaceOperator(const ThisType& )=delete;

  virtual inline void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  // dump system matrix into file
  void print(const std::string& filename="interface_matrix.dat") const
  {
    std::ofstream ofs(filename);
    const auto rows(op_.matrix().rows());
    auto count(decltype(rows){0});
    for(auto row=decltype(rows){0};row!=rows;++row)
    {
      while(count<(op_.matrix().numNonZeros()*(row+1)))
      {
        const auto entry(op_.matrix().realValue(count));
        const auto value(entry.first);
        const auto col(entry.second);
        if((std::abs(value)>1.e-13)&&(col>-1))
          ofs<<row+1<<" "<<col+1<<" "<<value<<std::endl;
        ++count;
      }
    }
  }

  inline const DiscreteSpaceType& domainSpace() const
  {
    return space_;
  }

  inline const DiscreteSpaceType& rangeSpace() const
  {
    return space_;
  }

  inline const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  template<typename BulkInterfaceGridMapperType,typename TimeProviderType>
  void assemble(const BulkInterfaceGridMapperType& mapper,const TimeProviderType& timeProvider) const
  {
    // allocate matrix
    DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    op_.reserve(stencil);
    op_.clear();
    // allocate local basis
    const auto blockSize(DiscreteSpaceType::localBlockSize);
    typedef typename DiscreteFunctionType::LocalFunctionType::RangeType LocalFunctionRangeType;
    std::vector<LocalFunctionRangeType> phi(space_.blockMapper().maxNumDofs()*blockSize );
    typedef typename DiscreteFunctionType::LocalFunctionType::JacobianRangeType LocalFunctionJacobianRangeType;
    std::vector<LocalFunctionJacobianRangeType> gradphi(space_.blockMapper().maxNumDofs()*blockSize);
    // extract dimensions
    constexpr auto worlddim(DiscreteSpaceType::GridType::dimensionworld);
    constexpr auto griddim(DiscreteSpaceType::GridType::dimension);
    constexpr auto rangedim(DiscreteSpaceType::FunctionSpaceType::dimRange);
    // define normal
    typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
    typedef typename DiscreteSpaceType::GridType::ctype ctype;
    typedef Normal<ctype,worlddim> NormalType;
    NormalType normal;
    typename NormalType::NormalVectorType normalVector;
    // assemble global matrix
    for(const auto entity:space_)
    {
      // compute normal
      const auto& faceLocalIdx(mapper.faceLocalIdxInterface2Bulk(space_.grid().leafIndexSet().index(entity)));
      normal(entity,normalVector,faceLocalIdx);
      // extract local matrix and basis functions
      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& baseSet(localMatrix.domainBasisFunctionSet());
      // assemble local \vec{A_m} (position)
      CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto qp:quadrature)
      {
        // evaluate basis functions and weight
        baseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        // fill \vec{A_m}
        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(std::size_t i=worlddim;i!=rowLocalSize;++i)
        {
          for(std::size_t j=worlddim;j!=columnLocalSize;++j)
          {
            RangeFieldType value(0.0);
            for(auto k=1;k!=rangedim;++k)
              value+=gradphi[i][k]*gradphi[j][k];
            value*=weight;
            localMatrix.add(i,j,value);
          }
        }
      }
      // assemble local \vec{N_m} (curvature_j-position_i)
      #if USE_LAGRANGE_QUADRATURE_POINTS
      auto pointSet(space_.template subDiscreteFunctionSpace<0>().lagrangePointSet(entity));
      const auto ptWeight(ReferenceElements<ctype,griddim>::general(entity.type()).volume()/static_cast<RangeFieldType>(pointSet.nop()));
      #else
      CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> pointSet(entity,2*space_.order()+1);
      #endif
      for(auto pt=0;pt!=pointSet.nop();++pt)
      {
        // evaluate basis functions and weight
        baseSet.evaluateAll(pointSet.point(pt),phi);
        #if USE_LAGRANGE_QUADRATURE_POINTS
        const auto weight(entity.geometry().integrationElement(pointSet.point(pt))*ptWeight);
        #else
        const auto weight(entity.geometry().integrationElement(pointSet.point(pt))*pointSet.weight(pt));
        #endif
        // fill \vec{N_m}
        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(std::size_t i=worlddim;i!=rowLocalSize;++i)
        {
          for(auto j=0;j!=columnLocalSize;++j)
          {
            RangeFieldType value(0.0);
            for(auto index=0;index!=worlddim;++index)
              value+=phi[i][index+1]*normalVector[index];
            value*=weight*phi[j][0];
            localMatrix.add(i,j,value);
            localMatrix.add(j,i,-1.0*value/(timeProvider.deltaT()));
          }
        }
      }
    }
  }

  private:
  const DiscreteSpaceType& space_;
  mutable LinearOperatorType op_;
};

}
}
#endif // DUNE_FEM_INTERFACEOPERATOR_HH
