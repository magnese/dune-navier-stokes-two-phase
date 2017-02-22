#ifndef DUNE_NORMAL_HH
#define DUNE_NORMAL_HH

namespace Dune
{

template<typename EntityType,int worlddim>
struct Normal
{};

// normal dim 2
template<typename EntityType>
struct Normal<EntityType,2>
{
  static auto compute(const EntityType& entity,unsigned int faceIdx)
  {
    typename EntityType::Geometry::GlobalCoordinate normalVector;
    normalVector[0]=(entity.geometry().corner(0)-entity.geometry().corner(1))[1];
    normalVector[1]=(entity.geometry().corner(1)-entity.geometry().corner(0))[0];
    normalVector/=normalVector.two_norm();
    if(faceIdx%2==1)
      normalVector*=(-1.0);
    return normalVector;
  }
};

// normal dim 3
template<typename EntityType>
struct Normal<EntityType,3>
{
  static auto compute(const EntityType& entity,unsigned int faceIdx)
  {
    typename EntityType::Geometry::GlobalCoordinate normalVector;
    const auto edge1(entity.geometry().corner(1)-entity.geometry().corner(0));
    const auto edge2(entity.geometry().corner(2)-entity.geometry().corner(0));
    normalVector[0]=edge1[1]*edge2[2]-edge1[2]*edge2[1];
    normalVector[1]=edge1[2]*edge2[0]-edge1[0]*edge2[2];
    normalVector[2]=edge1[0]*edge2[1]-edge1[1]*edge2[0];
    normalVector/=normalVector.two_norm();
    if(faceIdx%2==1)
      normalVector*=(-1.0);
    return normalVector;
  }
};

// helper function
template<typename EntityType>
inline auto computeNormal(const EntityType& entity,unsigned int faceIdx)
{
  return Normal<EntityType,EntityType::Geometry::coorddimension>::compute(entity,faceIdx);
}

}

#endif // DUNE_NORMAL_HH
