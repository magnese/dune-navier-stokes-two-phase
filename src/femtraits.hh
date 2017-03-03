#ifndef DUNE_FEM_FEMTRAITS_HH
#define DUNE_FEM_FEMTRAITS_HH

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/discontinuousgalerkin/lagrange.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/tuplediscretefunction.hh>

#include "extendedtuplediscretefunction.hh"

namespace Dune
{
namespace Fem
{

template<typename CoupledMeshManagerImp>
struct FemTraits
{
  // extract grids and grid parts
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef typename CoupledMeshManagerType::BulkGridType BulkGridType;
  typedef typename CoupledMeshManagerType::BulkGridPartType BulkGridPartType;
  typedef typename CoupledMeshManagerType::BulkInnerGridPartType BulkInnerGridPartType;
  typedef typename CoupledMeshManagerType::BulkOuterGridPartType BulkOuterGridPartType;
  typedef typename CoupledMeshManagerType::InterfaceGridType InterfaceGridType;
  typedef typename CoupledMeshManagerType::InterfaceGridPartType InterfaceGridPartType;
  // define continuos spaces
  typedef FunctionSpace<double,double,BulkGridType::dimensionworld,BulkGridType::dimensionworld> VelocityContinuosSpaceType;
  typedef FunctionSpace<double,double,BulkGridType::dimensionworld,1> PressureContinuosSpaceType;
  typedef FunctionSpace<double,double,BulkGridType::dimensionworld,BulkGridType::dimensionworld> BulkDisplacementContinuosSpaceType;
  typedef FunctionSpace<double,double,BulkGridType::dimensionworld,1> PhysicalCoefficientContinuosSpaceType;
  typedef FunctionSpace<double,double,InterfaceGridType::dimensionworld,InterfaceGridType::dimensionworld> DisplacementContinuosSpaceType;
  typedef FunctionSpace<double,double,InterfaceGridType::dimensionworld,1> CurvatureContinuosSpaceType;
  // define discrete spaces
  typedef LagrangeDiscreteFunctionSpace<VelocityContinuosSpaceType,BulkGridPartType,2> VelocityDiscreteSpaceType;
  #if PRESSURE_SPACE_TYPE == 0
  typedef LagrangeDiscontinuousGalerkinSpace<PressureContinuosSpaceType,BulkGridPartType,0> Pressure0DiscreteSpaceType;
  typedef Pressure0DiscreteSpaceType PressureDiscreteSpaceType;
  #elif PRESSURE_SPACE_TYPE == 1
  typedef LagrangeDiscreteFunctionSpace<PressureContinuosSpaceType,BulkGridPartType,1> Pressure0DiscreteSpaceType;
  typedef Pressure0DiscreteSpaceType PressureDiscreteSpaceType;
  #elif PRESSURE_SPACE_TYPE == 2
  typedef LagrangeDiscreteFunctionSpace<PressureContinuosSpaceType,BulkGridPartType,1> Pressure0DiscreteSpaceType;
  typedef LagrangeDiscontinuousGalerkinSpace<PressureContinuosSpaceType,BulkGridPartType,0> Pressure1DiscreteSpaceType;
  typedef LagrangeDiscontinuousGalerkinSpace<PressureContinuosSpaceType,BulkGridPartType,1> PressureDiscreteSpaceType;
  #endif
  typedef LagrangeDiscreteFunctionSpace<BulkDisplacementContinuosSpaceType,BulkGridPartType,1> BulkDisplacementDiscreteSpaceType;
  typedef LagrangeDiscontinuousGalerkinSpace<PhysicalCoefficientContinuosSpaceType,BulkGridPartType,0> PhysicalCoefficientDiscreteSpaceType;
  typedef LagrangeDiscreteFunctionSpace<DisplacementContinuosSpaceType,InterfaceGridPartType,1> DisplacementDiscreteSpaceType;
  typedef LagrangeDiscreteFunctionSpace<CurvatureContinuosSpaceType,InterfaceGridPartType,1> CurvatureDiscreteSpaceType;
  // define discrete functions
  typedef AdaptiveDiscreteFunction<VelocityDiscreteSpaceType> VelocityDiscreteFunctionType;
  typedef AdaptiveDiscreteFunction<Pressure0DiscreteSpaceType> Pressure0DiscreteFunctionType;
  #if PRESSURE_SPACE_TYPE == 0 || PRESSURE_SPACE_TYPE == 1
  typedef ExtendedTupleDiscreteFunction<VelocityDiscreteFunctionType,Pressure0DiscreteFunctionType> BulkDiscreteFunctionType;
  #elif PRESSURE_SPACE_TYPE == 2
  typedef AdaptiveDiscreteFunction<Pressure1DiscreteSpaceType> Pressure1DiscreteFunctionType;
  typedef ExtendedTupleDiscreteFunction<VelocityDiscreteFunctionType,Pressure0DiscreteFunctionType,Pressure1DiscreteFunctionType>
    BulkDiscreteFunctionType;
  #endif
  typedef AdaptiveDiscreteFunction<PressureDiscreteSpaceType> PressureDiscreteFunctionType;
  typedef AdaptiveDiscreteFunction<BulkDisplacementDiscreteSpaceType> BulkDisplacementDiscreteFunctionType;
  typedef AdaptiveDiscreteFunction<PhysicalCoefficientDiscreteSpaceType> PhysicalCoefficientDiscreteFunctionType;
  typedef AdaptiveDiscreteFunction<CurvatureDiscreteSpaceType> CurvatureDiscreteFunctionType;
  typedef AdaptiveDiscreteFunction<DisplacementDiscreteSpaceType> DisplacementDiscreteFunctionType;
  typedef TupleDiscreteFunction<CurvatureDiscreteFunctionType,DisplacementDiscreteFunctionType> InterfaceDiscreteFunctionType;
};

}
}

#endif // DUNE_FEM_FEMTRAITS_HH
