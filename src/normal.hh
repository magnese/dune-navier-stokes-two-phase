#ifndef DUNE_NORMAL_HH
#define DUNE_NORMAL_HH

#include <dune/common/fvector.hh>

namespace Dune
{

template<class ctype,int worlddim>
class Normal
{};

// normal dim 2
template<class ctype>
class Normal<ctype,2>
{
  public:
  typedef FieldVector<ctype,2> NormalVectorType;

  // compute the normal vector of a given entity
  template<class EntityType>
  void operator()(const EntityType& entity,NormalVectorType& normalVector,const unsigned int& faceIdx)
  {
    normalVector[0]=(entity.geometry().corner(0)-entity.geometry().corner(1))[1];
    normalVector[1]=(entity.geometry().corner(1)-entity.geometry().corner(0))[0];
    normalVector/=normalVector.two_norm();
    if(faceIdx%2==1)
      normalVector*=(-1.0);
  }
};

// normal dim 3
template<class ctype>
class Normal<ctype,3>
{
  public:
  typedef FieldVector<ctype,3> NormalVectorType;

  // compute the normal vector of a given entity
  template<class EntityType>
  void operator()(const EntityType& entity,NormalVectorType& normalVector,const unsigned int& faceIdx)
  {
    NormalVectorType edge1(entity.geometry().corner(1)-entity.geometry().corner(0));
    NormalVectorType edge2(entity.geometry().corner(2)-entity.geometry().corner(0));
    normalVector[0]=edge1[1]*edge2[2]-edge1[2]*edge2[1];
    normalVector[1]=edge1[2]*edge2[0]-edge1[0]*edge2[2];
    normalVector[2]=edge1[0]*edge2[1]-edge1[1]*edge2[0];
    normalVector/=normalVector.two_norm();
    if(faceIdx%2==1)
      normalVector*=(-1.0);
  }
};

}

#endif // DUNE_NORMAL_HH
