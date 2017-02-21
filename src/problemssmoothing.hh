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
    bc_.addAllBoundaryIDs();
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
