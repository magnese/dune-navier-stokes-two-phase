#ifndef DUNE_FEM_VERTEXFUNCTION_HH
#define DUNE_FEM_VERTEXFUNCTION_HH

// C++ includes
#include <string>
#include <algorithm>

// dune includes
#include <dune/common/exceptions.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/fem/gridpart/leafgridpart.hh>

#include <dune/grid/geometrygrid/coordfunction.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
//#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>

namespace Dune
{
namespace Fem
{

template<class GridImp>
class VertexFunction:public DiscreteCoordFunction<typename GridImp::ctype,GridImp::dimensionworld,VertexFunction<GridImp>>
{
  public:
  typedef GridImp GridType;
  typedef VertexFunction<GridType> ThisType;
  typedef typename GridType::ctype ctype;
  static constexpr auto griddim=GridType::dimension;
  static constexpr auto worlddim=GridType::dimensionworld;
  typedef DiscreteCoordFunction<ctype,worlddim,ThisType> BaseType;

  typedef typename BaseType::RangeVector RangeVectorType;

  typedef LeafGridPart<GridType> GridPartType;

  typedef FunctionSpace<ctype,ctype,worlddim,worlddim> ContinuousSpaceType;
  typedef LagrangeDiscreteFunctionSpace<ContinuousSpaceType,GridPartType,1,CachingStorage> DiscreteSpaceType;
  //typedef AdaptiveDiscreteFunction<DiscreteSpaceType> DiscreteFunctionType;
  typedef ISTLBlockVectorDiscreteFunction<DiscreteSpaceType> DiscreteFunctionType;

  typedef typename GridType::template Codim<0>::Entity HostEntityType;
  typedef typename GridType::template Codim<griddim>::Entity HostVertexType;

  explicit VertexFunction(GridType& grid):
    gridpart_(grid),space_(gridpart_),coord_("coordinate function",space_)
  {
    initialize(grid);
  }

  VertexFunction(const ThisType& )=delete;

  inline GridPartType& gridPart() const
  {
    return gridpart_;
  }
  inline DiscreteSpaceType& space() const
  {
    return space_;
  }
  inline DiscreteFunctionType& discreteFunction() const
  {
    return coord_;
  }

  void initialize(GridType& grid) const
  {
    // fill the vertices coordinates with the grid vertices
    const auto localBlockSize(DiscreteSpaceType::localBlockSize);
    for(const auto entity:entities(discreteFunction()))
    {
      auto localCoord(discreteFunction().localFunction(entity));
      const auto numLocalBlocks(entity.geometry().corners());
      std::size_t row(0);
      for(auto localIdx=0;localIdx!=numLocalBlocks;++localIdx)
      {
        const auto x(entity.geometry().corner(localIdx));
        for(auto l=0;l!=localBlockSize;++l,++row)
          localCoord[row]=x[l];
      }
    }
  }

  const ThisType& operator=(const std::vector<double>& vtx) const
  {
    std::copy(vtx.begin(),vtx.end(),discreteFunction().dbegin());
    return *this;
  }
  template<class DFT>
  const ThisType& operator=(const DFT& vtx) const
  {
    std::copy(vtx.dbegin(),vtx.dend(),discreteFunction().dbegin());
    return *this;
  }

  const ThisType& operator+=(const std::vector<double>& vtx) const
  {
    auto it(vtx.begin());
    for(auto& dof:dofs(discreteFunction()))
    {
      dof+=*it;
      ++it;
    }
    return *this;
  }
  template<class DFT>
  const ThisType& operator+=(const DFT& vtx) const
  {
    auto it(vtx.dbegin());
    for(auto& dof:dofs(discreteFunction()))
    {
      dof+=*it;
      ++it;
    }
    return *this;
  }

  const ThisType& operator-=(const std::vector<double>& vtx) const
  {
    auto it(vtx.begin());
    for(auto& dof:dofs(discreteFunction()))
    {
      dof-=*it;
      ++it;
    }
    return *this;
  }
  template<class DFT>
  const ThisType& operator-=(const DFT& vtx) const
  {
    auto it(vtx.dbegin());
    for(auto& dof:dofs(discreteFunction()))
    {
      dof-=*it;
      ++it;
    }
    return *this;
  }

  inline void evaluate(const HostEntityType& entity,unsigned int corner,RangeVectorType& y) const
  {
    const auto& referenceElement(ReferenceElements<ctype,griddim>::general((gridPart().template begin<0>())->type()));
    discreteFunction().localFunction(entity).evaluate(referenceElement.position(corner,griddim),y);
  }

  inline void evaluate(const HostVertexType& vertex,unsigned int corner,RangeVectorType& y) const
  {
    coord_.evaluate(vertex.geometry().corner(corner),y);;
  }

  inline void adapt() const
  {
    DUNE_THROW(NotImplemented,"VertexFunction::adapt() not implemented");
  }

  private:
  mutable GridPartType gridpart_;
  const DiscreteSpaceType space_;
  mutable DiscreteFunctionType coord_;
};

}
}
#endif // DUNE_FEM_VERTEXFUNCTION_HH
