#ifndef DUNE_FEM_NULLRHS_HH
#define DUNE_FEM_NULLRHS_HH

#include <fstream>
#include <string>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionImp,typename... Args>
class NullRHS
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef NullRHS<DiscreteFunctionType,Args...> ThisType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;

  template<typename... Argss>
  explicit NullRHS(const DiscreteSpaceType& space,const Argss&... ):
    space_(space),rhs_("null RHS",space_)
  {}

  NullRHS(const ThisType& )=delete;

  DiscreteFunctionType& rhs() const
  {
    return rhs_;
  }

  // dump rhs vector into file
  void print(const std::string& filename="null_rhs.dat") const
  {
    std::ofstream ofs(filename);
    rhs_.print(ofs);
  }

  double norm() const
  {
    return 0;
  }

  template<typename... Argss>
  void assemble(const Argss&... ) const
  {
    rhs_.clear();
  }

  private:
  const DiscreteSpaceType& space_;
  mutable DiscreteFunctionType rhs_;
};

}
}

#endif // DUNE_FEM_NULLRHS_HH
