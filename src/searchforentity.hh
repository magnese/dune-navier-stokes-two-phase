#ifndef DUNE_FEM_SEARCHFORENTITY_HH
#define DUNE_FEM_SEARCHFORENTITY_HH

#include <limits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

namespace Dune
{
namespace Fem
{

template<typename EntityType,typename GlobalCoordinateType>
FieldVector<double,EntityType::dimension+1> barycentricCoordinates(const EntityType& entity,const GlobalCoordinateType& x)
{
  FieldVector<double,EntityType::dimension+1> coor(1);
  const auto xLocal(entity.geometry().local(x));
  for(auto i=decltype(x.size()){0};i!=x.size();++i)
  {
    coor[0]-=xLocal[i];
    coor[i+1]=xLocal[i];
  }
  return coor;
}

template<typename GridPartType,typename EntityType,typename GlobalCoordinateType>
unsigned int searchForEntity(const GridPartType& gridPart,EntityType&& entity,const GlobalCoordinateType& x)
{
  constexpr double toll(1.e-12);
  // starting form entity, find the entity which contains x
  unsigned int numIterations(0);
  bool found(false);
  while(!found)
  {
    found=true;
    // compute barycentric coordinates
    const auto coor(barycentricCoordinates(entity,x));
    // check if x is inside the entity
    double minCoor(std::numeric_limits<double>::max());
    std::size_t faceIdx(0);
    for(auto i=decltype(coor.size()){0};i!=coor.size();++i)
      if(coor[i]<-toll)
        if(coor[i]<minCoor)
        {
          minCoor=coor[i];
          faceIdx=EntityType::dimension-i;
          found=false;
        }
    // search the point x in the neighbouring entity closer to x
    if(!found)
      for(const auto& intersection:intersections(gridPart,entity))
        if(static_cast<std::size_t>(intersection.indexInInside())==faceIdx)
        {
          if(intersection.boundary())
            DUNE_THROW(InvalidStateException,"ERROR: searchForEntity is trying to go outside the domain --> domain not convex!");
          entity=intersection.outside();
          break;
        }
    ++numIterations;
  }
  return numIterations;
}

}
}

#endif // DUNE_FEM_SEARCHFORENTITY_HH
