#ifndef DUNE_FEM_PROBLEMSSMOOTHING_HH
#define DUNE_FEM_PROBLEMSSMOOTHING_HH

#include "boundarycondition.hh"

namespace Dune
{
namespace Fem
{

template<typename DiscreteSpaceImp,typename CoupledMeshManagerImp>
class ParallelepipedGeometry
{
  public:
  typedef DiscreteSpaceImp DiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef ParallelepipedGeometry<DiscreteSpaceType,CoupledMeshManagerType> ThisType;

  // define BC
  typedef FreeSlipCondition<DiscreteSpaceType,DiscreteSpaceType,CoupledMeshManagerType> BoundaryConditionType;

  explicit ParallelepipedGeometry(CoupledMeshManagerType& meshManager):
    bc_(meshManager)
  {
    typedef typename BoundaryConditionType::DomainType DomainType;
    typedef typename BoundaryConditionType::EntityType EntityType;
    typedef typename BoundaryConditionType::RangeType RangeType;

    bc_.addBC(2,[](const DomainType& ,double ,const EntityType& )
                  {
                    RangeType value(1.0);
                    value[1]=0.0;
                    return value;
                  });
    bc_.addBC(3,[](const DomainType& ,double ,const EntityType& )
                  {
                    RangeType value(1.0);
                    value[0]=0.0;
                    return value;
                  });
    bc_.addBC(4,[](const DomainType& ,double ,const EntityType& )
                  {
                    RangeType value(1.0);
                    value[1]=0.0;
                    return value;
                  });
    bc_.addBC(5,[](const DomainType& ,double ,const EntityType& )
                  {
                    RangeType value(1.0);
                    value[0]=0.0;
                    return value;
                  });
    bc_.addBC(6,[](const DomainType& ,double ,const EntityType& )
                  {
                    RangeType value(1.0);
                    value[EntityType::Geometry::coorddimension-1]=0.0;
                    return value;
                  });
    bc_.addBC(7,[](const DomainType& ,double ,const EntityType& )
                  {
                    RangeType value(1.0);
                    value[EntityType::Geometry::coorddimension-1]=0.0;
                    return value;
                  });
  }

  BoundaryConditionType& bc()
  {
    return bc_;
  }

  private:
  BoundaryConditionType bc_;
};

}
}

#endif // DUNE_FEM_PROBLEMSSMOOTHING_HH
