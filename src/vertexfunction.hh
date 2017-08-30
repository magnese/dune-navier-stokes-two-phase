#ifndef DUNE_FEM_VERTEXFUNCTION_HH
#define DUNE_FEM_VERTEXFUNCTION_HH

#include <dune/common/exceptions.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/geometrygrid/coordfunction.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>

namespace Dune
{
namespace Fem
{

template<typename GridImp>
class VertexFunction:public DiscreteCoordFunction<typename GridImp::ctype,GridImp::dimensionworld,VertexFunction<GridImp>>
{
  public:
  typedef GridImp GridType;
  typedef VertexFunction<GridType> ThisType;
  typedef typename GridType::ctype ctype;
  static constexpr unsigned int griddim=GridType::dimension;
  static constexpr unsigned int worlddim=GridType::dimensionworld;
  typedef DiscreteCoordFunction<ctype,worlddim,ThisType> BaseType;
  typedef typename BaseType::RangeVector RangeVectorType;

  typedef LeafGridPart<GridType> GridPartType;

  typedef FunctionSpace<ctype,ctype,worlddim,worlddim> ContinuousSpaceType;
  typedef LagrangeDiscreteFunctionSpace<ContinuousSpaceType,GridPartType,1> DiscreteSpaceType;
  typedef AdaptiveDiscreteFunction<DiscreteSpaceType> DiscreteFunctionType;

  typedef typename GridType::template Codim<0>::Entity HostEntityType;
  typedef typename GridType::template Codim<griddim>::Entity HostVertexType;

  explicit VertexFunction(GridType& grid):
    gridpart_(grid),space_(gridpart_),coord_("coordinates",space_),localcoord_(coord_)
  {
    initialize(grid);
  }

  VertexFunction(const ThisType& )=delete;

  ThisType& operator=(const ThisType& other)
  {
    operator=(other.discreteFunction());
    return *this;
  }

  DiscreteFunctionType& discreteFunction()
  {
    return coord_;
  }
  const DiscreteFunctionType& discreteFunction() const
  {
    return coord_;
  }

  void initialize(const GridType& grid)
  {
    // fill the vertices coordinates with the grid vertices
    LocalContribution<DiscreteFunctionType,Assembly::Set> localCoord(coord_);
    constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
    for(const auto& entity:entities(coord_))
    {
      localCoord.bind(entity);
      const auto numLocalBlocks(entity.geometry().corners());
      std::size_t row(0);
      for(auto localIdx=decltype(numLocalBlocks){0};localIdx!=numLocalBlocks;++localIdx)
      {
        const auto x(entity.geometry().corner(localIdx));
        for(auto l=decltype(localBlockSize){0};l!=localBlockSize;++l,++row)
          localCoord[row]=x[l];
      }
      localCoord.unbind();
    }
  }

  template<typename DF>
  ThisType& operator=(const DF& df)
  {
    coord_.assign(df);
    return *this;
  }

  template<typename DF>
  ThisType& operator+=(const DF& df)
  {
    coord_+=df;
    return *this;
  }

  template<typename DF>
  ThisType& operator-=(const DF& df)
  {
    coord_-=df;
    return *this;
  }

  void evaluate(const HostEntityType& entity,unsigned int corner,RangeVectorType& y) const
  {
    const auto refElement(referenceElement(entity.geometry()));
    localcoord_.init(entity);
    localcoord_.evaluate(refElement.position(corner,griddim),y);
  }

  void evaluate(const HostVertexType& vertex,unsigned int ,RangeVectorType& y) const
  {
    coord_.evaluate(vertex.geometry().center(),y);
  }

  void adapt()
  {
    DUNE_THROW(NotImplemented,"VertexFunction::adapt() not implemented");
  }

  private:
  GridPartType gridpart_;
  const DiscreteSpaceType space_;
  DiscreteFunctionType coord_;
  mutable ConstLocalDiscreteFunction<DiscreteFunctionType> localcoord_;
};

}
}
#endif // DUNE_FEM_VERTEXFUNCTION_HH
