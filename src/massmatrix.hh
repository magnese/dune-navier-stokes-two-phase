#ifndef DUNE_FEM_MASSMATRIX_HH
#define DUNE_FEM_MASSMATRIX_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <string>
#include <fstream>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DomainFunctionImp, typename RangeFunctionImp,
         template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class MassMatrix:public Operator<DomainFunctionImp,RangeFunctionImp>
{
  public:
  typedef DomainFunctionImp DomainFunctionType;
  typedef RangeFunctionImp RangeFunctionType;
  typedef LinearOperatorImp<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef MassMatrix<DomainFunctionType,RangeFunctionType,LinearOperatorImp> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef DomainFunctionType DiscreteFunctionType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef DomainSpaceType DiscreteSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  explicit MassMatrix(const DiscreteSpaceType& space):
    space_(space),op_("mass matrix",space_,space_)
  {}

  MassMatrix(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="mass_matrix.dat") const
  {
    std::ofstream ofs(filename);
    op_.matrix().print(ofs);
  }

  const DomainSpaceType& domainSpace() const
  {
    return space_;
  }

  const RangeSpaceType& rangeSpace() const
  {
    return space_;
  }

  void assemble()
  {
    DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    op_.reserve(stencil);
    op_.clear();

    constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
    typedef typename DiscreteFunctionType::LocalFunctionType::RangeType LocalFunctionRangeType;
    std::vector<LocalFunctionRangeType> phi(space_.blockMapper().maxNumDofs()*localBlockSize);

    // perform a grid walkthrough and assemble the global matrix
    for(const auto& entity:space_)
    {
      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& baseSet(localMatrix.domainBasisFunctionSet());
      CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto& qp:quadrature)
      {
        baseSet.evaluateAll(qp,phi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        const auto localSize(localMatrix.rows());
        for(auto i=decltype(localSize){0};i!=localSize;++i)
        {
          for(auto j=decltype(localSize){0};j!=localSize;++j)
          {
            const auto value((phi[i]*phi[j])*weight);
            localMatrix.add(i,j,value);
          }
        }
      }
    }

  }

  private:
  const DiscreteSpaceType& space_;
  LinearOperatorType op_;
};

}
}

#endif // DUNE_FEM_MASSMATRIX_HH
