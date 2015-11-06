#ifndef DUNE_FEM_COMPUTE_HH
#define DUNE_FEM_COMPUTE_HH

// C++ includes
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

// dune includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>

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

  //create time provider
  FixedStepTimeProvider<> timeProvider;

  // get parameters
  const double endTime(Dune::Fem::Parameter::getValue<double>("EndTime",0.1)+0.1*timeProvider.deltaT());
  const bool computeErrors(Parameter::getValue<bool>("ComputeErrors",0)&&femScheme.problem().hasExactSolution());

  // rebuild all quantities if the mesh is changed
  fluidState.update();

  // get initial conditions of the problem (if it is NOT time dependent the two solutions are garbage)
  femScheme.problem().velocityIC(fluidState.velocity());
  femScheme.problem().pressureIC(fluidState.pressureDump());

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

  // solve
  for(;timeProvider.time()<=endTime;timeProvider.next())
  {
    // print time
    std::cout<<std::endl<<"Time step "<<timeProvider.timeStep()<<" (time = "<<timeProvider.time()<<" s)."<<std::endl;

    // do one step
    femScheme(timeProvider);

    // dump solution
    fluidState.dumpBulkSolutions(timeProvider);
    fluidState.dumpInterfaceSolutions(timeProvider);

    // compute bulk errors
    if(computeErrors)
    {
      // interpolate exact solutions
      typename FluidStateType::VelocityDiscreteFunctionType velocityExactSolution("velocity exact solution",fluidState.velocitySpace());
      femScheme.problem().velocitySolution(velocityExactSolution,timeProvider.time());
      typename FluidStateType::PressureDumpDiscreteFunctionType pressureExactSolution("pressure exact solution",
                                                                                      fluidState.pressureDumpSpace());
      femScheme.problem().pressureSolution(pressureExactSolution,timeProvider.time());
      // compute L2 errors
      L2Norm<typename FluidStateType::BulkGridPartType> norm(fluidState.bulkGridPart());
      errors[1]+=pow(norm.distance(velocityExactSolution,fluidState.velocity()),2);
      errors[3]+=pow(norm.distance(pressureExactSolution,fluidState.pressureDump()),2);
      // compute Linfinity errors
      auto exactVelocityIt(velocityExactSolution.dbegin());
      for(auto velocityIt=fluidState.velocity().dbegin();velocityIt!=fluidState.velocity().dend();++velocityIt,++exactVelocityIt)
        errors[2]=std::max(errors[2],std::abs(*velocityIt-*exactVelocityIt));
      auto exactPressureIt(pressureExactSolution.dbegin());
      for(auto pressureIt=fluidState.pressureDump().dbegin();pressureIt!=fluidState.pressureDump().dend();++pressureIt,++exactPressureIt)
        errors[4]=std::max(errors[4],std::abs(*pressureIt-*exactPressureIt));
    }

    // update interface grid
    fluidState.interfaceGrid().coordFunction()+=fluidState.displacement();

    // compute radius error and L2 pressure error not interpolated
    #if PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 3
    if(computeErrors)
    {
      // compute radius error
      const auto interfaceGridLeafView(fluidState.interfaceGrid().leafGridView());
      for(const auto vertex:vertices(interfaceGridLeafView))
        errors[0]=std::max(errors[0],std::abs(femScheme.problem().exactRadius(timeProvider.time())-vertex.geometry().center().two_norm()));
      // store an inner and an outer entity, needed for indicator function
      const auto innerEntity(fluidState.bulkGrid().entity(fluidState.meshManager().mapper().entitySeedInterface2Bulk(0)));
      auto intersectionIt(fluidState.bulkGridPart().ibegin(innerEntity));
      while(intersectionIt->indexInInside()!=fluidState.meshManager().mapper().faceLocalIdxInterface2Bulk(0))
        ++intersectionIt;
      const auto outerEntity(intersectionIt->outside());
      // compute L2 pressure error not interpolated
      for(const auto entity:fluidState.pressureDumpSpace())
      {
        auto localPressureDump(fluidState.pressureDump().localFunction(entity));
        CachingQuadrature<typename FluidStateType::BulkGridPartType,0> pointSet(entity,13);
        for(auto pt=0;pt!=pointSet.nop();++pt)
        {
          const auto localPoint(pointSet.point(pt));
          typename FluidStateType::PressureDumpDiscreteFunctionType::RangeType valueDiscretePressure;
          localPressureDump.evaluate(localPoint,valueDiscretePressure);
          const auto globalPoint(entity.geometry().global(localPoint));
          typename FluidStateType::PressureDiscreteFunctionType::RangeType valueExactPressure;
          if(globalPoint.two_norm()>femScheme.problem().exactRadius(timeProvider.time()))
            valueExactPressure=femScheme.problem().pressureSolution()(globalPoint,timeProvider.time(),outerEntity);
          else
            valueExactPressure=femScheme.problem().pressureSolution()(globalPoint,timeProvider.time(),innerEntity);
          const auto weight(entity.geometry().integrationElement(localPoint)*pointSet.weight(pt));
          errors[5]+=pow(valueDiscretePressure-valueExactPressure,2)*weight;
        }
      }
    }
    #endif

    // keep also the original fluid state in order to interpolate the velocity onto the new grid
    auto oldFluidState(fluidState);

    // preform mesh smoothing
    fluidState.bulkGrid().coordFunction()+=meshSmoothing(fluidState.displacement());

    // perform remesh (if needed)
    fluidState.meshManager().remesh();

    // interpolate velocity onto the new grid (if necessary)
    if(!(femScheme.problem().isDensityNull()))
    {
      // rebuild all quantities if the mesh is changed
      fluidState.update();
      // interpolate velocity onto the new grid
      const auto velocityLocalBlockSize(FluidStateType::VelocityDiscreteSpaceType::localBlockSize);
      const auto& velocitySpace(fluidState.velocitySpace());
      for(const auto entity:velocitySpace)
      {
        auto localVelocity(fluidState.velocity().localFunction(entity));
        const auto& lagrangePointSet(velocitySpace.lagrangePointSet(entity));
        std::size_t row(0);
        for(auto pt=0;pt!=lagrangePointSet.nop();++pt)
        {
          typename FluidStateType::VelocityDiscreteFunctionType::RangeType temp;
          oldFluidState.velocity().evaluate(entity.geometry().global(lagrangePointSet.point(pt)),temp);
          for(auto l=0;l!=velocityLocalBlockSize;++l,++row)
            localVelocity[row]=temp[l];
        }
      }
    }
  }

  // print errors
  if(computeErrors)
  {
    errors[1]=pow(errors[1]*timeProvider.deltaT(),0.5);
    errors[3]=pow(errors[3]*timeProvider.deltaT(),0.5);
    errors[5]=pow(errors[5]*timeProvider.deltaT(),0.5);
    std::cout<<std::endl<<"Errors:"<<std::endl;
    std::cout.precision(5);
    std::cout<<std::scientific;
    #if PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 3
    std::cout<<"||X-x||_{L^oo} = "<<errors[0]<<std::endl;
    #endif
    std::cout<<"||U-Iu||_{L^2} = "<<errors[1]<<std::endl;
    std::cout<<"||U-Iu||_{L^oo} = "<<errors[2]<<std::endl;
    std::cout<<"||P-Ip||_{L^2} = "<<errors[3]<<std::endl;
    std::cout<<"||P-Ip||_{L^oo} = "<<errors[4]<<std::endl;
    #if PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 3
    std::cout<<"||P-p||_{L^2} = "<<errors[5]<<std::endl;
    #endif
    std::cout.unsetf(std::ios_base::floatfield);
  }
}

}
}

#endif // DUNE_FEM_COMPUTE_HH
