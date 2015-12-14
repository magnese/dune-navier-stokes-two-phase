#ifndef DUNE_FEM_INTERFACEDISPLACEMENTRHS_HH
#define DUNE_FEM_INTERFACEDISPLACEMENTRHS_HH

// deprecation warning
#warning ("WARNING : interfacedisplacementrhs.hh is deprecated")

#include <fstream>
#include <string>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionImp>
class InterfaceDisplacementRHS
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef InterfaceDisplacementRHS<DiscreteFunctionType> ThisType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;

  explicit InterfaceDisplacementRHS(const DiscreteSpaceType& space):
    space_(space),rhs_("interface displacement RHS",space_)
  {}

  InterfaceDisplacementRHS(const ThisType& )=delete;

  DiscreteFunctionType& rhs() const
  {
    return rhs_;
  }

  // dump rhs vector into file
  void print(const std::string& filename="interface_displacement_rhs.dat") const
  {
    std::ofstream ofs(filename);
    rhs_.print(ofs);
  }

  double norm() const
  {
    return sqrt(rhs_.scalarProductDofs(rhs_));
  }

  template<typename OperatorType>
  void assemble(const OperatorType& op) const
  {
    rhs_.clear();

    const auto& coord(space_.grid().coordFunction().discreteFunction());
    op(coord,rhs_);

    for(auto& dof:dofs(rhs_))
      dof*=-1.0;
  }

  private:
  const DiscreteSpaceType& space_;
  mutable DiscreteFunctionType rhs_;
};

}
}

#endif // DUNE_FEM_INTERFACEDISPLACEMENTRHS_HH
