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
  static constexpr int velocitySpacePolynomialOrder=2;
  typedef LagrangeDiscreteFunctionSpace<VelocityContinuosSpaceType,BulkGridPartType,velocitySpacePolynomialOrder> VelocityDiscreteSpaceType;
  #if PRESSURE_SPACE_TYPE == 0
  static constexpr int pressureSpacePolynomialOrder=0;
  typedef LagrangeDiscontinuousGalerkinSpace<PressureContinuosSpaceType,BulkGridPartType,pressureSpacePolynomialOrder> Pressure0DiscreteSpaceType;
  typedef Pressure0DiscreteSpaceType PressureDiscreteSpaceType;
  #elif PRESSURE_SPACE_TYPE == 1
  static constexpr int pressureSpacePolynomialOrder=1;
  typedef LagrangeDiscreteFunctionSpace<PressureContinuosSpaceType,BulkGridPartType,pressureSpacePolynomialOrder> Pressure0DiscreteSpaceType;
  typedef Pressure0DiscreteSpaceType PressureDiscreteSpaceType;
  #elif PRESSURE_SPACE_TYPE == 2
  typedef LagrangeDiscreteFunctionSpace<PressureContinuosSpaceType,BulkGridPartType,1> Pressure0DiscreteSpaceType;
  typedef LagrangeDiscontinuousGalerkinSpace<PressureContinuosSpaceType,BulkGridPartType,0> Pressure1DiscreteSpaceType;
  static constexpr int pressureSpacePolynomialOrder=1;
  typedef LagrangeDiscontinuousGalerkinSpace<PressureContinuosSpaceType,BulkGridPartType,pressureSpacePolynomialOrder> PressureDiscreteSpaceType;
  #elif PRESSURE_SPACE_TYPE == 3
  typedef LagrangeDiscreteFunctionSpace<PressureContinuosSpaceType,BulkInnerGridPartType,1> Pressure0DiscreteSpaceType;
  typedef LagrangeDiscreteFunctionSpace<PressureContinuosSpaceType,BulkOuterGridPartType,1> Pressure1DiscreteSpaceType;
  static constexpr int pressureSpacePolynomialOrder=1;
  typedef LagrangeDiscontinuousGalerkinSpace<PressureContinuosSpaceType,BulkGridPartType,pressureSpacePolynomialOrder> PressureDiscreteSpaceType;
  #endif
  typedef LagrangeDiscreteFunctionSpace<BulkDisplacementContinuosSpaceType,BulkGridPartType,1> BulkDisplacementDiscreteSpaceType;
  static constexpr int physicalCoefficientSpacePolynomialOrder=0;
  typedef LagrangeDiscontinuousGalerkinSpace<PhysicalCoefficientContinuosSpaceType,BulkGridPartType,physicalCoefficientSpacePolynomialOrder> PhysicalCoefficientDiscreteSpaceType;
  typedef LagrangeDiscreteFunctionSpace<DisplacementContinuosSpaceType,InterfaceGridPartType,1> DisplacementDiscreteSpaceType;
  typedef LagrangeDiscreteFunctionSpace<CurvatureContinuosSpaceType,InterfaceGridPartType,1> CurvatureDiscreteSpaceType;
  // define discrete functions
  typedef AdaptiveDiscreteFunction<VelocityDiscreteSpaceType> VelocityDiscreteFunctionType;
  typedef AdaptiveDiscreteFunction<Pressure0DiscreteSpaceType> Pressure0DiscreteFunctionType;
  #if !USE_EXTENDED_PRESSURE_SPACE
  typedef ExtendedTupleDiscreteFunction<VelocityDiscreteFunctionType,Pressure0DiscreteFunctionType> BulkDiscreteFunctionType;
  #else
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
