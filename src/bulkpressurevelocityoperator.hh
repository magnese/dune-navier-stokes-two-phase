#ifndef DUNE_FEM_BULKPRESSUREVELOCITYOPERATOR_HH
#define DUNE_FEM_BULKPRESSUREVELOCITYOPERATOR_HH

#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <fstream>
#include <string>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DomainFunctionImp,typename RangeFunctionImp,
         template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class BulkPressureVelocityOperator:public Operator<DomainFunctionImp,RangeFunctionImp>
{
  public:
  typedef DomainFunctionImp DomainFunctionType;
  typedef RangeFunctionImp RangeFunctionType;
  typedef LinearOperatorImp<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef BulkPressureVelocityOperator<DomainFunctionImp,RangeFunctionImp,LinearOperatorImp> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  explicit BulkPressureVelocityOperator(const DomainSpaceType& domainSpace,const RangeSpaceType& rangeSpace):
    domainspace_(domainSpace),rangespace_(rangeSpace),op_("bulk pressure velocity operator",domainspace_,rangespace_)
  {}

  BulkPressureVelocityOperator(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="pressure_velocity_matrix.dat") const
  {
    std::ofstream ofs(Parameter::getValue<std::string>("fem.prefix",".")+"/"+filename);
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

    // clear matrix
    op_.clear();

    // allocate basis vectors
    constexpr std::size_t domainLocalBlockSize(DomainSpaceType::localBlockSize);
    constexpr std::size_t rangeLocalBlockSize(RangeSpaceType::localBlockSize);
    typedef typename DomainFunctionType::LocalFunctionType::RangeType DomainRangeType;
    std::vector<DomainRangeType> phi(domainspace_.blockMapper().maxNumDofs()*domainLocalBlockSize);
    typedef typename RangeFunctionType::LocalFunctionType::JacobianRangeType RangeJacobianRangeType;
    std::vector<RangeJacobianRangeType> gradphi(rangespace_.blockMapper().maxNumDofs()*rangeLocalBlockSize);

    // perform a grid walkthrough and assemble the global matrix
    for(const auto& entity:domainspace_)
    {
      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& domainBaseSet(localMatrix.domainBasisFunctionSet());
      const auto& rangeBaseSet(localMatrix.rangeBasisFunctionSet());
      const CachingQuadrature<typename DomainSpaceType::GridPartType,0> quadrature(entity,2*domainspace_.order()+1);
      for(const auto& qp:quadrature)
      {
        domainBaseSet.evaluateAll(qp,phi);
        rangeBaseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

        const auto rowLocalSize(localMatrix.rows());
        const auto columnLocalSize(localMatrix.columns());
        for(auto i=decltype(rowLocalSize){0};i!=rowLocalSize;++i)
        {
          for(auto j=decltype(columnLocalSize){0};j!=columnLocalSize;++j)
          {
            typename RangeSpaceType::RangeFieldType value(0.0);
            for(auto k=decltype(rangeLocalBlockSize)(0);k!=rangeLocalBlockSize;++k)
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
  LinearOperatorType op_;
};

}
}

#endif // DUNE_FEM_BULKPRESSUREVELOCITYOPERATOR_HH
