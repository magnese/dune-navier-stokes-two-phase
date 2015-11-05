#ifndef DUNE_FEM_FEM_TRAITS_HH
#define DUNE_FEM_FEM_TRAITS_HH

// dune includes
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#if PRESSURE_SPACE_TYPE != 1
#include <dune/fem/space/discontinuousgalerkin/lagrange.hh>
#endif
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/function/blockvectorfunction.hh>

namespace Dune
{
namespace Fem
{

template<typename BulkGridImp,typename InterfaceGridImp>
struct FemTraits
{
  // define grids
  typedef BulkGridImp BulkGridType;
  typedef LeafGridPart<BulkGridType> BulkGridPartType;
  typedef InterfaceGridImp InterfaceGridType;
  typedef LeafGridPart<InterfaceGridType> InterfaceGridPartType;
  // define continuos spaces
  typedef FunctionSpace<double,double,BulkGridType::dimensionworld,BulkGridType::dimensionworld> VelocityContinuosSpaceType;
  typedef FunctionSpace<double,double,BulkGridType::dimensionworld,1> PressureContinuosSpaceType;
  typedef FunctionSpace<double,double,InterfaceGridType::dimensionworld,InterfaceGridType::dimensionworld> DisplacementContinuosSpaceType;
  typedef FunctionSpace<double,double,InterfaceGridType::dimensionworld,1> CurvatureContinuosSpaceType;
  // define discrete spaces
  typedef LagrangeDiscreteFunctionSpace<VelocityContinuosSpaceType,BulkGridPartType,2> VelocityDiscreteSpaceType;
  #if PRESSURE_SPACE_TYPE == 0
  typedef LagrangeDiscontinuousGalerkinSpace<PressureContinuosSpaceType,BulkGridPartType,0> PressureDiscreteSpaceType;
  typedef PressureDiscreteSpaceType PressureDumpDiscreteSpaceType;
  #elif PRESSURE_SPACE_TYPE == 1
  typedef LagrangeDiscreteFunctionSpace<PressureContinuosSpaceType,BulkGridPartType,1> PressureDiscreteSpaceType;
  typedef PressureDiscreteSpaceType PressureDumpDiscreteSpaceType;
  #else
  typedef LagrangeDiscreteFunctionSpace<PressureContinuosSpaceType,BulkGridPartType,1> PressureDiscreteSpaceType;
  typedef LagrangeDiscontinuousGalerkinSpace<PressureContinuosSpaceType,BulkGridPartType,0> PressureAdditionalDiscreteSpaceType;
  typedef LagrangeDiscontinuousGalerkinSpace<PressureContinuosSpaceType,BulkGridPartType,1> PressureDumpDiscreteSpaceType;
  #endif
  typedef LagrangeDiscreteFunctionSpace<DisplacementContinuosSpaceType,InterfaceGridPartType,1> DisplacementDiscreteSpaceType;
  typedef LagrangeDiscreteFunctionSpace<CurvatureContinuosSpaceType,InterfaceGridPartType,1> CurvatureDiscreteSpaceType;
  #if PRESSURE_SPACE_TYPE !=2
  typedef TupleDiscreteFunctionSpace<VelocityDiscreteSpaceType,PressureDiscreteSpaceType> BulkDiscreteSpaceType;
  #else
  typedef TupleDiscreteFunctionSpace<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,PressureAdditionalDiscreteSpaceType>
    BulkDiscreteSpaceType;
  #endif
  typedef TupleDiscreteFunctionSpace<CurvatureDiscreteSpaceType,DisplacementDiscreteSpaceType> InterfaceDiscreteSpaceType;
  // define discrete functions
  typedef ISTLBlockVectorDiscreteFunction<VelocityDiscreteSpaceType> VelocityDiscreteFunctionType;
  typedef ISTLBlockVectorDiscreteFunction<PressureDiscreteSpaceType> PressureDiscreteFunctionType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef ISTLBlockVectorDiscreteFunction<PressureAdditionalDiscreteSpaceType> PressureAdditionalDiscreteFunctionType;
  #endif
  typedef ISTLBlockVectorDiscreteFunction<PressureDumpDiscreteSpaceType> PressureDumpDiscreteFunctionType;
  typedef ISTLBlockVectorDiscreteFunction<CurvatureDiscreteSpaceType> CurvatureDiscreteFunctionType;
  typedef ISTLBlockVectorDiscreteFunction<DisplacementDiscreteSpaceType> DisplacementDiscreteFunctionType;
  typedef ISTLBlockVectorDiscreteFunction<BulkDiscreteSpaceType> BulkDiscreteFunctionType;
  typedef ISTLBlockVectorDiscreteFunction<InterfaceDiscreteSpaceType> InterfaceDiscreteFunctionType;
};

}
}

#endif
