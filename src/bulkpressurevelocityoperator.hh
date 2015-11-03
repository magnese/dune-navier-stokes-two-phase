#ifndef DUNE_FEM_BULKPRESSUREVELOCITYOPERATOR_HH
#define DUNE_FEM_BULKPRESSUREVELOCITYOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <string>
#include <fstream>

namespace Dune
{
namespace Fem
{

template<class LinearOperatorImp>
class BulkPressureVelocityOperator:public Operator<typename LinearOperatorImp::DomainFunctionType,
                                                   typename LinearOperatorImp::RangeFunctionType>
{
  public:
  typedef LinearOperatorImp LinearOperatorType;
  typedef typename LinearOperatorType::DomainFunctionType DomainFunctionType;
  typedef typename LinearOperatorType::RangeFunctionType RangeFunctionType;
  typedef BulkPressureVelocityOperator<LinearOperatorType> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  explicit BulkPressureVelocityOperator(const DomainSpaceType& domainSpace,const RangeSpaceType& rangeSpace):
    domainspace_(domainSpace),rangespace_(rangeSpace),op_("bulk pressure velocity operator",domainspace_,rangespace_)
  {}

  BulkPressureVelocityOperator(const ThisType& )=delete;

  inline LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  // apply the operator
  virtual inline void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  // dump system matrix into file
  inline void print(const std::string& filename="pressure_velocity_matrix.dat") const
  {
    std::ofstream ofs(filename);
    op_.matrix().print(ofs,1);
    ofs.close();
  }

  inline const DomainSpaceType& domainSpace() const
  {
    return domainspace_;
  }

  inline const RangeSpaceType& rangeSpace() const
  {
    return rangespace_;
  }

  void assemble() const
  {
    DiagonalAndNeighborStencil<DomainSpaceType,RangeSpaceType> stencil(domainspace_,rangespace_);
    op_.reserve(stencil);
    op_.clear();

    const auto domainLocalBlockSize(DomainSpaceType::localBlockSize);
    const auto rangeLocalBlockSize(RangeSpaceType::localBlockSize);
    typedef typename DomainFunctionType::LocalFunctionType::RangeType DomainRangeType;
    std::vector<DomainRangeType> phi(domainspace_.blockMapper().maxNumDofs()*domainLocalBlockSize);
    typedef typename RangeFunctionType::LocalFunctionType::JacobianRangeType RangeJacobianRangeType;
    std::vector<RangeJacobianRangeType> gradphi(rangespace_.blockMapper().maxNumDofs()*rangeLocalBlockSize);

    // perform a grid walkthrough and assemble the global matrix
    for(const auto entity:domainspace_)
    {
      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& domainBaseSet(localMatrix.domainBasisFunctionSet());
      const auto& rangeBaseSet(localMatrix.rangeBasisFunctionSet());
      CachingQuadrature<typename DomainSpaceType::GridPartType,0> pointSet(entity,2*domainspace_.order()+1);
      for(auto pt=0;pt!=pointSet.nop();++pt)
      {
        // evaluate the jacobians of all basis functions
        domainBaseSet.evaluateAll(pointSet.point(pt),phi);
        rangeBaseSet.jacobianAll(pointSet.point(pt),gradphi);
        const auto weight(entity.geometry().integrationElement(pointSet.point(pt))*pointSet.weight(pt));

        const auto rowLocalSize(localMatrix.rows());
        const auto columnLocalSize(localMatrix.columns());
        for(auto i=0;i!=rowLocalSize;++i)
        {
          for(auto j=0;j!=columnLocalSize;++j)
          {
            typename RangeSpaceType::RangeFieldType value(0.0);
            for(auto k=0;k!=rangeLocalBlockSize;++k)
              value+=gradphi[i][k][k];
            value*=(-1.0*weight*phi[j][0]);
            localMatrix.add(i,j,value);
          }
        }
      }
    }

  }

  private:
  const DomainSpaceType& domainspace_;
  const RangeSpaceType& rangespace_;
  mutable LinearOperatorType op_;
};

}
}

#endif // DUNE_FEM_BULKPRESSUREVELOCITYOPERATOR_HH
