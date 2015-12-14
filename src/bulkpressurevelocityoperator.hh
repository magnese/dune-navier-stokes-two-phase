#ifndef DUNE_FEM_BULKPRESSUREVELOCITYOPERATOR_HH
#define DUNE_FEM_BULKPRESSUREVELOCITYOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <string>
#include <fstream>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename LinearOperatorImp>
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

  LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  // apply the operator
  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  // dump system matrix into file
  void print(const std::string& filename="pressure_velocity_matrix.dat") const
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

  const DomainSpaceType& domainSpace() const
  {
    return domainspace_;
  }

  const RangeSpaceType& rangeSpace() const
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
      CachingQuadrature<typename DomainSpaceType::GridPartType,0> quadrature(entity,2*domainspace_.order()+1);
      for(const auto qp:quadrature)
      {
        // evaluate the jacobians of all basis functions
        domainBaseSet.evaluateAll(qp,phi);
        rangeBaseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

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
