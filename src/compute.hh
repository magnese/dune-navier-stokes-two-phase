#ifndef DUNE_FEM_COMPUTE_HH
#define DUNE_FEM_COMPUTE_HH

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>

#include "differentmeshlocalevaluator.hh"
#include "sortedview.hh"

#if PROBLEM_NUMBER == 7
#include "bubblestatistics.hh"
#endif

namespace Dune
{
namespace Fem
{

template<typename FemSchemeType,typename MeshSmoothingType>
void compute(FemSchemeType& femScheme,MeshSmoothingType& meshSmoothing,std::vector<double>& errors)
{
  // define fluid state
  typedef typename FemSchemeType::FluidStateType FluidStateType;
  auto& fluidState(femScheme.fluidState());

  // create time provider
  FixedStepTimeProvider<> timeProvider;

  // get parameters
  const double endTime(Parameter::getValue<double>("EndTime",0.1)+0.1*timeProvider.deltaT());
  const bool computeErrors(Parameter::getValue<bool>("ComputeErrors",0)&&femScheme.problem().hasExactSolution());

  // rebuild all quantities if the mesh is changed
  fluidState.update();

  // get initial conditions of the problem
  femScheme.problem().velocityIC(fluidState.velocity());
  femScheme.problem().pressureIC(fluidState.pressureDump());
  #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
  if(!femScheme.problem().isDensityNull())
    femScheme.problem().rho(fluidState.rho());
  #endif

  // dump bulk solution at t0 and advance time provider
  if(femScheme.problem().isTimeDependent())
  {
    // compute initial curvature of the interface
    femScheme.computeInterface(timeProvider);
    // dump solution on file
    fluidState.dumpBulkSolutions(timeProvider);
    fluidState.dumpInterfaceSolutions(timeProvider);
    // advance time provider
    timeProvider.next();
  }

  // clear errors
  errors.clear();
  errors.resize(6,0.0);

  // compute bulk bounding box which is needed by interpolation 2
  #if INTERPOLATION_TYPE == 2
  typename FluidStateType::CoupledMeshManagerType::BulkBoundingBoxType bulkBoundingBox;
  if(!(femScheme.problem().isDensityNull()))
  {
    bulkBoundingBox=fluidState.meshManager().bulkBoundingBox();
    std::cout<<"Bulk bounding box : ["<<bulkBoundingBox.first<<"] ; ["<<bulkBoundingBox.second<<"]\n";
  }
  #endif

  #if PROBLEM_NUMBER == 7
  BubbleStatistics bubbleStatistics;
  bubbleStatistics.add(fluidState,timeProvider);
  #endif

  // solve
  for(;timeProvider.time()<=endTime;timeProvider.next())
  {
    // print time
    std::cout<<"\nTime step "<<timeProvider.timeStep()<<" (time = "<<timeProvider.time()<<" s).\n";
    Timer timerStep(false);
    timerStep.start();

    // do one step
    femScheme(timeProvider,meshSmoothing);

    // dump solution
    fluidState.dumpBulkSolutions(timeProvider);
    fluidState.dumpInterfaceSolutions(timeProvider);

    #if PROBLEM_NUMBER == 7
    bubbleStatistics.add(fluidState,timeProvider);
    #endif

    // compute bulk errors
    if(computeErrors)
    {
      // interpolate exact solutions
      typename FluidStateType::VelocityDiscreteFunctionType velocityExactSolution("velocity exact solution",fluidState.velocitySpace());
      femScheme.problem().velocitySolution(velocityExactSolution,timeProvider.time());
      typename FluidStateType::PressureDumpDiscreteFunctionType pressureExactSolution("pressure exact solution",
                                                                                      fluidState.pressureDumpSpace());
      femScheme.problem().pressureSolution(pressureExactSolution,timeProvider.time());
      // compute L2 interpolated velocity error
      L2Norm<typename FluidStateType::BulkGridPartType> normL2(fluidState.bulkGridPart());
      errors[1]+=std::pow(normL2.distance(velocityExactSolution,fluidState.velocity()),2);
      // compute H1 interpolated velocity error
      H1Norm<typename FluidStateType::BulkGridPartType> normH1(fluidState.bulkGridPart());
      errors[2]+=std::pow(normH1.distance(velocityExactSolution,fluidState.velocity()),2);
      // compute Linfinity interpolated velocity error
      auto exactVelocityIt(velocityExactSolution.dbegin());
      for(const auto& dof:dofs(fluidState.velocity()))
        errors[3]=std::max(errors[3],std::abs(dof-*(exactVelocityIt++)));
      #if PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 3 || PROBLEM_NUMBER == 8 || PROBLEM_NUMBER == 9 || PROBLEM_NUMBER == 10
      // store an inner and an outer entity, needed for indicator function
      const auto interfaceEntity(*(fluidState.interfaceGridPart().template begin<0>());
      const auto intersection(fluidState.meshManager().correspondingInnerBulkIntersection(interfaceEntity));
      const auto innerEntity(intersection.inside());
      const auto outerEntity(intersection.outside());
      // compute L2 pressure error
      for(const auto& entity:fluidState.pressureDumpSpace())
      {
        auto localPressureDump(fluidState.pressureDump().localFunction(entity));
        constexpr unsigned int order(FluidStateType::BulkGridType::dimensionworld<3?13:10);
        CachingQuadrature<typename FluidStateType::BulkGridPartType,0> quadrature(entity,order);
        for(const auto& qp:quadrature)
        {
          const auto localPoint(qp.position());
          typename FluidStateType::PressureDumpDiscreteFunctionType::RangeType value;
          localPressureDump.evaluate(localPoint,value);
          const auto globalPoint(entity.geometry().global(localPoint));
          if(globalPoint.two_norm()>femScheme.problem().exactRadius(timeProvider.time()))
            value-=femScheme.problem().pressureSolution().function()(globalPoint,timeProvider.time(),outerEntity);
          else
            value-=femScheme.problem().pressureSolution().function()(globalPoint,timeProvider.time(),innerEntity);
          const auto weight(entity.geometry().integrationElement(localPoint)*qp.weight());
          errors[4]+=value.two_norm2()*weight;
        }
      }
      // compute Linfinity radius error
      fluidState.interfaceGrid().coordFunction()+=fluidState.displacement();
      for(const auto& vertex:vertices(fluidState.interfaceGridPart()))
        errors[0]=std::max(errors[0],std::abs(femScheme.problem().exactRadius(timeProvider.time())-vertex.geometry().center().two_norm()));
      fluidState.interfaceGrid().coordFunction()-=fluidState.displacement();
      #endif
      // compute Linfinity interpolated pressure error
      auto exactPressureIt(pressureExactSolution.dbegin());
      for(const auto& dof:dofs(fluidState.pressureDump()))
        errors[5]=std::max(errors[5],std::abs(dof-*(exactPressureIt++)));
    }

    // update bulk and interface grid
    const double nonLinearSolverType(Parameter::getValue<int>("NonLinearSolverType",0));
    const bool useALE(nonLinearSolverType>1);
    fluidState.interfaceGrid().coordFunction()+=fluidState.displacement();
    if(femScheme.problem().isDensityNull()||(!useALE))
      meshSmoothing.computeBulkDisplacement();
    fluidState.bulkGrid().coordFunction()+=fluidState.bulkDisplacement();

    // perform remesh and keep also the original fluid state to interpolate the velocity onto the new grid
    auto oldFluidState(fluidState);
    const auto remeshPerformed(fluidState.meshManager().remesh());
    auto interpolationNeeded((!femScheme.problem().isDensityNull())&&remeshPerformed);

    // check if the smoothing has modified the bulk mesh (only when velocity interpolation is needed)
    if((!femScheme.problem().isDensityNull())&&(!remeshPerformed)&&(!useALE))
    {
      const double nullTolerance(Parameter::getValue<double>("NullTolerance",1.e-12));
      // check if the bulk displacement is not null
      for(const auto& dof:dofs(fluidState.bulkDisplacement()))
        if(std::abs(dof)>nullTolerance)
        {
          interpolationNeeded=true;
          break;
        }
    }

    // interpolate velocity and rho onto the new grid
    if(interpolationNeeded)
    {
      Timer timerInterpolation(false);
      timerInterpolation.start();
      // rebuild all quantities if the mesh is changed
      fluidState.update();
      // set old bulk grid and old velocity to the correct values
      if(remeshPerformed)
        oldFluidState.bulkGrid().coordFunction()-=oldFluidState.bulkDisplacement();
      else
      {
        // deep copy of the old mesh manager from the new mesh manager to have independent fluid states
        oldFluidState.meshManager().deepCopy(fluidState.meshManager());
        oldFluidState.init();
        // set old bulk grid, old velocity and old rho to the correct values
        oldFluidState.bulkGrid().coordFunction()-=fluidState.bulkDisplacement();
        oldFluidState.velocity().assign(fluidState.velocity());
        #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
        oldFluidState.rho().assign(fluidState.rho());
        #endif
      }
      #if INTERPOLATION_TYPE == 0
      constexpr bool useBarycentricEntitySearch(false);
      auto& newBulkGridPart(fluidState.bulkGridPart());
      #elif INTERPOLATION_TYPE == 1
      constexpr bool useBarycentricEntitySearch(true);
      auto& newBulkGridPart(fluidState.bulkGridPart());
      #else
      constexpr bool useBarycentricEntitySearch(true);
      SortedView<typename FluidStateType::BulkGridType> newBulkGridPart(fluidState.bulkGrid(),bulkBoundingBox);
      #endif
      constexpr std::size_t velocityLocalBlockSize(FluidStateType::VelocityDiscreteSpaceType::localBlockSize);
      std::vector<typename FluidStateType::VelocityDiscreteFunctionType::RangeFieldType> velocityInterpolatedDOFs;
      velocityInterpolatedDOFs.reserve(fluidState.velocitySpace().blockMapper().maxNumDofs()*velocityLocalBlockSize);
      DifferentMeshLocalEvaluator<typename FluidStateType::VelocityDiscreteFunctionType,useBarycentricEntitySearch>
        localOldVelocity(oldFluidState.velocity());
      #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
      constexpr std::size_t rhoLocalBlockSize(FluidStateType::PhysicalCoefficientDiscreteSpaceType::localBlockSize);
      std::vector<typename FluidStateType::PhysicalCoefficientDiscreteFunctionType::RangeFieldType> rhoInterpolatedDOFs;
      rhoInterpolatedDOFs.reserve(fluidState.rhoSpace().blockMapper().maxNumDofs()*rhoLocalBlockSize);
      DifferentMeshLocalEvaluator<typename FluidStateType::PhysicalCoefficientDiscreteFunctionType,useBarycentricEntitySearch>
        localOldRho(oldFluidState.rho());
      #endif
      for(const auto& entity:elements(newBulkGridPart))
      {
        localOldVelocity.init(entity);
        velocityInterpolatedDOFs.resize(fluidState.velocitySpace().basisFunctionSet(entity).size());
        fluidState.velocitySpace().interpolation(entity)(localOldVelocity,velocityInterpolatedDOFs);
        fluidState.velocity().setLocalDofs(entity,velocityInterpolatedDOFs);
        #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
        localOldRho.init(entity,localOldVelocity.oldEntity());
        rhoInterpolatedDOFs.resize(fluidState.rhoSpace().basisFunctionSet(entity).size());
        fluidState.rhoSpace().interpolation(entity)(localOldRho,rhoInterpolatedDOFs);
        fluidState.rho().setLocalDofs(entity,rhoInterpolatedDOFs);
        #endif
      }
      timerInterpolation.stop();
      if(useBarycentricEntitySearch)
        std::cout<<"Average research depth : "<<
          static_cast<double>(localOldVelocity.searchIterations())/static_cast<double>(fluidState.velocitySpace().blockMapper().size())
          <<".\n";
      std::cout<<"Velocity interpolation time: "<<timerInterpolation.elapsed()<<" seconds.\n";
    }
    timerStep.stop();

    std::cout<<"Full timestep time: "<<timerStep.elapsed()<<" seconds.\n";
  }

  // print errors
  if(computeErrors)
  {
    errors[1]=std::pow(errors[1]*timeProvider.deltaT(),0.5);
    errors[2]=std::pow(errors[2]*timeProvider.deltaT(),0.5);
    errors[4]=std::pow(errors[4]*timeProvider.deltaT(),0.5);
    std::cout<<"\nErrors:\n";
    std::cout.precision(5);
    std::cout<<std::scientific;
    #if PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 3 || PROBLEM_NUMBER == 8 || PROBLEM_NUMBER == 9 || PROBLEM_NUMBER == 10
    std::cout<<"||X-x||_{L^oo} = "<<errors[0]<<"\n";
    #endif
    std::cout<<"||U-Iu||_{L^2} = "<<errors[1]<<"\n";
    std::cout<<"||U-Iu||_{H^1} = "<<errors[2]<<"\n";
    std::cout<<"||U-Iu||_{L^oo} = "<<errors[3]<<"\n";
    #if PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 3 || PROBLEM_NUMBER == 8 || PROBLEM_NUMBER == 9 || PROBLEM_NUMBER == 10
    std::cout<<"||P-p||_{L^2} = "<<errors[4]<<"\n";
    #endif
    std::cout<<"||P-Ip||_{L^oo} = "<<errors[5]<<"\n";
    std::cout.unsetf(std::ios_base::floatfield);
  }
}

}
}

#endif // DUNE_FEM_COMPUTE_HH
