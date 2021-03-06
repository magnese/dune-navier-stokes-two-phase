#ifndef DUNE_FEM_MASSMATRIX_HH
#define DUNE_FEM_MASSMATRIX_HH

#include <dune/fem/io/io.hh>
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

  void print(const std::string& filename="mass_matrix.dat",unsigned int offset=0) const
  {
    const std::string& path(Parameter::getValue<std::string>("fem.prefix","."));
    if(!directoryExists(path))
      createDirectory(path);
    std::ofstream ofs(path+"/"+filename);
    op_.exportMatrix().print(ofs,offset);
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
    // allocate matrix
    Stencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    for(const auto& entity:space_)
      stencil.fill(entity,entity);
    op_.reserve(stencil);

    // clear matrix
    op_.clear();

    // allocate basis vector
    std::vector<typename DiscreteFunctionType::RangeType> phi(space_.maxNumDofs());

    // perform a grid walkthrough and assemble the global matrix
    for(const auto& entity:space_)
    {
      auto localMatrix(op_.localMatrix(entity,entity));
      const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto& qp:quadrature)
      {
        localMatrix.domainBasisFunctionSet().evaluateAll(qp,phi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        const auto localSize(localMatrix.rows());
        for(auto i=decltype(localSize){0};i!=localSize;++i)
          for(auto j=decltype(localSize){0};j!=localSize;++j)
            localMatrix.add(i,j,(phi[i]*phi[j])*weight);
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
