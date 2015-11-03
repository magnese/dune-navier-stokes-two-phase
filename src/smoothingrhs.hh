#ifndef DUNE_FEM_SMOOTHINGRHS_HH
#define DUNE_FEM_SMOOTHINGRHS_HH

#include <string>
#include <fstream>

namespace Dune
{
namespace Fem
{

template<class DiscreteFunctionImp>
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

  inline DiscreteFunctionType& rhs() const
  {
    return rhs_;
  }

  // dump rhs vector into file
  inline void print(const std::string& filename="smoothing_rhs.dat") const
  {
    std::ofstream ofs(filename);
    rhs_.print(ofs);
    ofs.close();
  }

  inline double norm() const
  {
    return sqrt(rhs_.scalarProductDofs(rhs_));
  }

  template<class T,class M>
  inline void assemble(const T& interfaceDisplacement,const M& bulkInterfaceGridMapper) const
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
