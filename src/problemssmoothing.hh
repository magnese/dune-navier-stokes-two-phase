#ifndef DUNE_FEM_PROBLEMSSMOOTHING_HH
#define DUNE_FEM_PROBLEMSSMOOTHING_HH

#include <dune/fem/function/common/localfunctionadapter.hh>

#include "boundarycondition.hh"

namespace Dune
{
namespace Fem
{

template<typename FluidStateImp>
class ParallelepipedGeometry
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef ParallelepipedGeometry<FluidStateType> ThisType;

  // define BC
  typedef typename FluidStateType::BulkDisplacementDiscreteSpaceType DiscreteSpaceType;
  typedef typename FluidStateType::CoupledMeshManagerType CoupledMeshManagerType;
  typedef FreeSlipCondition<DiscreteSpaceType,CoupledMeshManagerType> BoundaryConditionType;

  typedef typename DiscreteSpaceType::DomainType DomainType;
  typedef typename DiscreteSpaceType::RangeType RangeType;
  typedef typename DiscreteSpaceType::EntityType EntityType;
  typedef LocalAnalyticalFunctionBinder<DiscreteSpaceType> LocalAnalyticalFunctionType;

  explicit ParallelepipedGeometry(FluidStateType& fluidState):
    bc_(fluidState.meshManager())
  {
    bc_.addBC(2,LocalAnalyticalFunctionType([](const DomainType& ,double ,const EntityType& )
                                            {
                                              RangeType value(1.0);
                                              value[1]=0.0;
                                              return value;
                                            }));
    bc_.addBC(3,LocalAnalyticalFunctionType([](const DomainType& ,double ,const EntityType& )
                                            {
                                              RangeType value(1.0);
                                              value[0]=0.0;
                                              return value;
                                            }));
    bc_.addBC(4,LocalAnalyticalFunctionType([](const DomainType& ,double ,const EntityType& )
                                            {
                                              RangeType value(1.0);
                                              value[1]=0.0;
                                              return value;
                                            }));
    bc_.addBC(5,LocalAnalyticalFunctionType([](const DomainType& ,double ,const EntityType& )
                                            {
                                              RangeType value(1.0);
                                              value[0]=0.0;
                                              return value;
                                            }));
    bc_.addBC(6,LocalAnalyticalFunctionType([](const DomainType& ,double ,const EntityType& )
                                            {
                                              RangeType value(1.0);
                                              value[EntityType::Geometry::coorddimension-1]=0.0;
                                              return value;
                                            }));
    bc_.addBC(7,LocalAnalyticalFunctionType([](const DomainType& ,double ,const EntityType& )
                                            {
                                              RangeType value(1.0);
                                              value[EntityType::Geometry::coorddimension-1]=0.0;
                                              return value;
                                            }));
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
