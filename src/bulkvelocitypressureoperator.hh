#ifndef DUNE_FEM_BULKVELOCITYPRESSUREOPERATOR_HH
#define DUNE_FEM_BULKVELOCITYPRESSUREOPERATOR_HH

#include <dune/fem/operator/linear/spoperator.hh>
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

template<typename DomainFunctionImp,typename RangeFunctionImp,
         template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class BulkVelocityPressureOperator:public Operator<DomainFunctionImp,RangeFunctionImp>
{
  public:
  typedef DomainFunctionImp DomainFunctionType;
  typedef RangeFunctionImp RangeFunctionType;
  typedef LinearOperatorImp<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef BulkVelocityPressureOperator<DomainFunctionType,RangeFunctionType,LinearOperatorImp> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  explicit BulkVelocityPressureOperator(const DomainSpaceType& domainSpace,const RangeSpaceType& rangeSpace):
    domainspace_(domainSpace),rangespace_(rangeSpace),op_("bulk velocity pressure operator",domainspace_,rangespace_)
  {}

  BulkVelocityPressureOperator(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="velocity_pressure_matrix.dat") const
  {
    std::ofstream ofs(filename);
    op_.matrix().print(ofs);
  }

  const DomainSpaceType& domainSpace() const
  {
    return domainspace_;
  }

  const RangeSpaceType& rangeSpace() const
  {
    return rangespace_;
  }

  void assemble()
  {
    DiagonalAndNeighborStencil<DomainSpaceType,RangeSpaceType> stencil(domainspace_,rangespace_);
    op_.reserve(stencil);
    op_.clear();

    const auto domainLocalBlockSize(DomainSpaceType::localBlockSize);
    const auto rangeLocalBlockSize(RangeSpaceType::localBlockSize);
    typedef typename DomainFunctionType::LocalFunctionType::JacobianRangeType DomainJacobianRangeType;
    std::vector<DomainJacobianRangeType> gradphi(domainspace_.blockMapper().maxNumDofs()*domainLocalBlockSize);
    typedef typename RangeFunctionType::LocalFunctionType::RangeType RangeRangeType;
    std::vector<RangeRangeType> phi(rangespace_.blockMapper().maxNumDofs()*rangeLocalBlockSize);

    for(const auto& entity:domainspace_)
    {
      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& domainBaseSet(localMatrix.domainBasisFunctionSet());
      const auto& rangeBaseSet(localMatrix.rangeBasisFunctionSet());

      CachingQuadrature<typename DomainSpaceType::GridPartType,0> quadrature(entity,2*domainspace_.order()+1);
      for(const auto& qp:quadrature)
      {
        rangeBaseSet.evaluateAll(qp,phi);
        domainBaseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(auto i=0;i!=rowLocalSize;++i)
        {
          for(auto j=0;j!=columnLocalSize;++j)
          {
            typename RangeSpaceType::RangeFieldType value(0.0);
            for(auto k=0;k!=domainLocalBlockSize;++k)
              value+=gradphi[j][k][k];
            value*=(-1.0*weight*phi[i][0]);
            localMatrix.add(i,j,value);
          }
        }
      }
    }

  }

  private:
  const DomainSpaceType& domainspace_;
  const RangeSpaceType& rangespace_;
  LinearOperatorType op_;
};

}
}

#endif // DUNE_FEM_BULKVELOCITYPRESSUREOPERATOR_HH
