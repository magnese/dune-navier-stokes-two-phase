#ifndef DUNE_FEM_INTERFACERHS_HH
#define DUNE_FEM_INTERFACERHS_HH

#include <fstream>
#include <string>
#include <algorithm>

namespace Dune
{
namespace Fem
{

template<class DiscreteFunctionImp>
class InterfaceRHS
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef InterfaceRHS<DiscreteFunctionType> ThisType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;

  explicit InterfaceRHS(const DiscreteSpaceType& space):
    space_(space),rhs_("interface RHS",space_)
  {}

  InterfaceRHS(const ThisType& )=delete;

  inline DiscreteFunctionType& rhs() const
  {
    return rhs_;
  }

  // dump rhs vector into file
  inline void print(const std::string& filename="interface_rhs.dat") const
  {
    std::ofstream ofs(filename);
    rhs_.print(ofs);
    ofs.close();
  }

  inline double norm() const
  {
    return sqrt(rhs_.scalarProductDofs(rhs_));
  }

  template<class OperatorType>
  void assemble(const OperatorType& op) const
  {
    rhs_.clear();

    DiscreteFunctionType temp("temp",space_);
    const auto space0Size(space_.template subDiscreteFunctionSpace<0>().size());
    auto lastIt(std::fill_n(temp.dbegin(),space0Size,0.0));

    const auto& coord(space_.grid().coordFunction().discreteFunction());
    std::copy(coord.dbegin(),coord.dend(),lastIt);

    op(temp,rhs_);

    std::fill_n(rhs_.dbegin(),space0Size,0.0);
  }

  private:
  const DiscreteSpaceType& space_;
  mutable DiscreteFunctionType rhs_;
};

}
}

#endif // DUNE_FEM_INTERFACERHS_HH
