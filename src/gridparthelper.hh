#ifndef DUNE_FEM_GRIDPARTHELPER_HH
#define DUNE_FEM_GRIDPARTHELPER_HH

#include <tuple>

#include <dune/fem/gridpart/filteredgridpart.hh>

namespace Dune
{
namespace Fem
{

template<typename GridPartType>
struct GridPartHelper
{
  typedef typename GridPartType::GridType::template Codim<0>::Entity EntityType;

  static bool contains(const GridPartType& ,const EntityType& )
  {
    return true;
  }
};

template<typename GridPartType,typename IndicatorFunctionType,bool UseFilteredIndexSet>
struct GridPartHelper<FilteredGridPart<GridPartType,IndicatorFunctionType,UseFilteredIndexSet>>
{
  typedef FilteredGridPart<GridPartType,IndicatorFunctionType,UseFilteredIndexSet> FilteredGridPartType;
  typedef typename FilteredGridPartType::GridType::template Codim<0>::Entity EntityType;

  static bool contains(const FilteredGridPartType& filteredGridPart,const EntityType& entity)
  {
    return filteredGridPart.contains(entity);
  }
};

template<typename EntityType,typename... GridPartsType>
auto isEntityContained(const EntityType& entity,const GridPartsType&... gridParts)
{
  return std::make_tuple(GridPartHelper<GridPartsType>::contains(gridParts,entity)...);
}

};
};

#endif // DUNE_FEM_GRIDPARTHELPER_HH
