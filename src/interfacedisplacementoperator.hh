#ifndef DUNE_FEM_INTERFACEDISPLACEMENTOPERATOR_HH
#define DUNE_FEM_INTERFACEDISPLACEMENTOPERATOR_HH

// deprecation warning
#warning ("WARNING : interfacedisplacementoperator.hh is deprecated")

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <vector>
#include <string>
#include <fstream>

namespace Dune
{
namespace Fem
{

template<typename LinearOperatorImp>
class InterfaceDisplacementOperator:public Operator<typename LinearOperatorImp::DomainFunctionType,
                                                    typename LinearOperatorImp::RangeFunctionType>
{
  public:
  typedef LinearOperatorImp LinearOperatorType;
  typedef typename LinearOperatorType::DomainFunctionType DomainFunctionType;
  typedef typename LinearOperatorType::RangeFunctionType RangeFunctionType;
  typedef DomainFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;
  typedef InterfaceDisplacementOperator<LinearOperatorType> ThisType;

  InterfaceDisplacementOperator(const ThisType& )=delete;

  // constructor
  explicit InterfaceDisplacementOperator(const DiscreteSpaceType& space):
    space_(space),op_("interface displacement operator",space_,space_)
  {}

  // apply the operator
  virtual inline void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  template<typename DT,typename RT>
  inline void operator()(const DT& u,RT& w) const
  {
    op_.apply(u,w);
  }

  // dump matrix into file
  void print(const std::string& filename="interface_displacement_matrix.dat") const
  {
    std::ofstream ofs(filename);
    op_.matrix().print(ofs,1);
    ofs.close();
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

  // perform a grid walkthrough and assemble the global matrix
  void assemble() const
  {
    DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    op_.reserve(stencil);
    op_.clear();

    const auto localBlockSize(DiscreteSpaceType::localBlockSize);
    typedef typename DiscreteFunctionType::LocalFunctionType::JacobianRangeType LocalFunctionJacobianRangeType;
    std::vector<LocalFunctionJacobianRangeType> gradphi(space_.blockMapper().maxNumDofs()*localBlockSize);

    constexpr auto rangedim(DiscreteSpaceType::FunctionSpaceType::dimRange);

    // perform a grid walkthrough and assemble the global matrix
    for(const auto entity:space_)
    {
      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& baseSet(localMatrix.domainBasisFunctionSet());

      CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> pointSet(entity,2*space_.order()+1);
      for(auto pt=0;pt!=pointSet.nop();++pt)
      {
        baseSet.jacobianAll(pointSet.point(pt),gradphi);
        const auto weight(entity.geometry().integrationElement(pointSet.point(pt))*pointSet.weight(pt));

        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(auto i=0;i!=rowLocalSize;++i)
        {
          for(auto j=0;j!=columnLocalSize;++j)
          {
            typename DiscreteSpaceType::RangeFieldType value(0.0);
            for(auto k=0;k!=rangedim;++k)
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
  mutable LinearOperatorType op_;
};

}
}
#endif // DUNE_FEM_INTERFACEDISPLACEMENTOPERATOR_HH
