#ifndef DUNE_FEM_INTERFACEDISPLACEMENTOPERATOR_HH
#define DUNE_FEM_INTERFACEDISPLACEMENTOPERATOR_HH

// deprecation warning
#warning ("WARNING : interfacedisplacementoperator.hh is deprecated")

#include <dune/fem/io/io.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/common/operator/linear/spoperator.hh>

#include <string>
#include <fstream>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionImp,template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class InterfaceDisplacementOperator:public Operator<DiscreteFunctionImp,DiscreteFunctionImp>
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef DiscreteFunctionType DomainFunctionType;
  typedef DiscreteFunctionType RangeFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  typedef LinearOperatorImp<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef typename LinearOperatorType::MatrixType MatrixType;
  typedef InterfaceDisplacementOperator<DiscreteFunctionType,LinearOperatorImp> ThisType;

  InterfaceDisplacementOperator(const ThisType& )=delete;

  // constructor
  explicit InterfaceDisplacementOperator(const DiscreteSpaceType& space):
    space_(space),op_("interface displacement operator",space_,space_)
  {}

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  template<typename DT,typename RT>
  void operator()(const DT& u,RT& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="interface_displacement_matrix.dat",unsigned int offset=0) const
  {
    const std::string& path(Parameter::getValue<std::string>("fem.prefix","."));
    if(!directoryExists(path))
      createDirectory(path);
    std::ofstream ofs(path+"/"+filename);
    op_.matrix().print(ofs,offset);
  }

  const DiscreteSpaceType& domainSpace() const
  {
    return space_;
  }

  const DiscreteSpaceType& rangeSpace() const
  {
    return space_;
  }

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  void assemble()
  {
    DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    op_.reserve(stencil);
    op_.clear();

    typedef typename DiscreteFunctionType::LocalFunctionType::JacobianRangeType LocalFunctionJacobianRangeType;
    std::vector<LocalFunctionJacobianRangeType> gradphi(space_.maxNumDofs());

    constexpr std::size_t rangedim(DiscreteSpaceType::FunctionSpaceType::dimRange);

    for(const auto& entity:space_)
    {
      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& baseSet(localMatrix.domainBasisFunctionSet());

      const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto& qp:quadrature)
      {
        baseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(auto i=decltype(rowLocalSize){0};i!=rowLocalSize;++i)
        {
          for(auto j=decltype(columnLocalSize){0};j!=columnLocalSize;++j)
          {
            typename DiscreteSpaceType::RangeFieldType value(0.0);
            for(auto k=decltype(rangedim){0};k!=rangedim;++k)
              value+=gradphi[i][k]*gradphi[j][k];
            value*=weight;
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
#endif // DUNE_FEM_INTERFACEDISPLACEMENTOPERATOR_HH
