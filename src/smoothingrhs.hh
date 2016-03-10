#ifndef DUNE_FEM_SMOOTHINGRHS_HH
#define DUNE_FEM_SMOOTHINGRHS_HH

#include <cmath>
#include <string>
#include <fstream>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionImp>
class SmoothingRHS
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef SmoothingRHS<DiscreteFunctionType> ThisType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;

  explicit SmoothingRHS(const DiscreteSpaceType& space):
    space_(space),rhs_("smoothing RHS",space_)
  {}

  SmoothingRHS(const ThisType& )=delete;

  DiscreteFunctionType& rhs() const
  {
    return rhs_;
  }

  // dump rhs vector into file
  void print(const std::string& filename="smoothing_rhs.dat") const
  {
    std::ofstream ofs(filename);
    rhs_.print(ofs);
  }

  double norm() const
  {
    return std::sqrt(rhs_.scalarProductDofs(rhs_));
  }

  template<typename T,typename M>
  void assemble(const T& interfaceDisplacement,const M& bulkInterfaceGridMapper) const
  {
    rhs_.clear();
    bulkInterfaceGridMapper.addInterfaceDF2BulkDF(interfaceDisplacement,rhs_);
  }

  private:
  const DiscreteSpaceType& space_;
  mutable DiscreteFunctionType rhs_;
};

}
}

#endif // DUNE_FEM_SMOOTHINGRHS_HH
