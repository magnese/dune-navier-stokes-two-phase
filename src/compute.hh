#ifndef DUNE_FEM_COMPUTE_HH
#define DUNE_FEM_COMPUTE_HH

// C++ includes
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

// dune includes
#include <dune/common/exceptions.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>

#include "searchforentity.hh"
#include "sortedview.hh"

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
    Timer timerStep(false);
    Timer timerInterpolation(false);

    // print time
    std::cout<<std::endl<<"Time step "<<timeProvider.timeStep()<<" (time = "<<timeProvider.time()<<" s)."<<std::endl;
    timerStep.start();

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
      // compute L2 interpolated velocity error
      L2Norm<typename FluidStateType::BulkGridPartType> normL2(fluidState.bulkGridPart());
      errors[1]+=std::pow(normL2.distance(velocityExactSolution,fluidState.velocity()),2);
      // compute H1 interpolated velocity error
      H1Norm<typename FluidStateType::BulkGridPartType> normH1(fluidState.bulkGridPart());
      errors[2]+=std::pow(normH1.distance(velocityExactSolution,fluidState.velocity()),2);
      // compute Linfinity interpolated velocity error
      auto exactVelocityIt(velocityExactSolution.dbegin());
      for(auto velocityIt=fluidState.velocity().dbegin();velocityIt!=fluidState.velocity().dend();++velocityIt,++exactVelocityIt)
        errors[3]=std::max(errors[3],std::abs(*velocityIt-*exactVelocityIt));
      // compute L2 pressure error
      #if PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 3
      // store an inner and an outer entity, needed for indicator function
      const auto innerEntity(fluidState.bulkGrid().entity(fluidState.meshManager().mapper().entitySeedInterface2Bulk(0)));
      auto intersectionIt(fluidState.bulkGridPart().ibegin(innerEntity));
      while(intersectionIt->indexInInside()!=fluidState.meshManager().mapper().faceLocalIdxInterface2Bulk(0))
        ++intersectionIt;
      const auto outerEntity(intersectionIt->outside());
      // compute error
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
            value-=femScheme.problem().pressureSolution()(globalPoint,timeProvider.time(),outerEntity);
          else
            value-=femScheme.problem().pressureSolution()(globalPoint,timeProvider.time(),innerEntity);
          const auto weight(entity.geometry().integrationElement(localPoint)*qp.weight());
          errors[4]+=value.two_norm2()*weight;
        }
      }
      #endif
      // compute Linfinity interpolated pressure error
      auto exactPressureIt(pressureExactSolution.dbegin());
      for(auto pressureIt=fluidState.pressureDump().dbegin();pressureIt!=fluidState.pressureDump().dend();++pressureIt,++exactPressureIt)
        errors[5]=std::max(errors[5],std::abs(*pressureIt-*exactPressureIt));
    }

    // update interface grid
    fluidState.interfaceGrid().coordFunction()+=fluidState.displacement();

    // compute Linfinity radius error
    if(computeErrors)
    {
      #if PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 3
      const auto interfaceGridLeafView(fluidState.interfaceGrid().leafGridView());
      for(const auto& vertex:vertices(interfaceGridLeafView))
        errors[0]=std::max(errors[0],std::abs(femScheme.problem().exactRadius(timeProvider.time())-vertex.geometry().center().two_norm()));
      #endif
    }

    // perform mesh smoothing
    if(meshSmoothing.isEnabled())
      meshSmoothing.apply();
    else
    {
      fluidState.bulkDisplacement().clear();
      fluidState.meshManager().mapper().addInterfaceDF2BulkDF(fluidState.displacement(),fluidState.bulkDisplacement());
    }
    fluidState.bulkGrid().coordFunction()+=fluidState.bulkDisplacement();

    // perform remesh and keep also the original fluid state to interpolate the velocity onto the new grid
    auto oldFluidState(fluidState);
    const auto remeshPerformed(fluidState.meshManager().remesh());
    timerInterpolation.start();

    // interpolate velocity onto the new grid
    if(!(femScheme.problem().isDensityNull()))
    {
      // rebuild all quantities if the mesh is changed
      fluidState.update();
      // restore old bulk grid
      //oldFluidState.interfaceGrid().coordFunction()-=oldFluidState.displacement();
      if(meshSmoothing.isEnabled()||!remeshPerformed)
        DUNE_THROW(InvalidStateException,"ERROR: for Navier-Stokes smoothing needs to be always OFF and remeshing always ON!");
      else
        oldFluidState.bulkGrid().coordFunction()-=oldFluidState.bulkDisplacement();
      // store if a dof has already been interpolated
      std::vector<bool> dofAlreadyInterpolated(fluidState.velocitySpace().blockMapper().size(),false);
      // interpolate velocity onto the new grid
      const auto velocityLocalBlockSize(FluidStateType::VelocityDiscreteSpaceType::localBlockSize);
      const auto& velocitySpace(fluidState.velocitySpace());
      #if INTERPOLATION_TYPE != 0
      auto oldEntity=(*(oldFluidState.velocitySpace().begin()));
      unsigned int averageResearchDepth(0);
      #endif
      #if INTERPOLATION_TYPE == 2
      SortedView<typename FluidStateType::BulkGridType> sortedView(fluidState.bulkGrid(),0.025,{0,0},{1,2});
      for(const auto& entity:elements(sortedView))
      #else
      for(const auto& entity:velocitySpace)
      #endif
      {
        auto localVelocity(fluidState.velocity().localFunction(entity));
        const auto& lagrangePointSet(velocitySpace.lagrangePointSet(entity));
        std::size_t row(0);
        std::vector<std::size_t> globalIdxs(velocitySpace.blockMapper().numDofs(entity));
        velocitySpace.blockMapper().map(entity,globalIdxs);
        for(auto pt=decltype(lagrangePointSet.nop()){0};pt!=lagrangePointSet.nop();++pt)
        {
          if(dofAlreadyInterpolated[globalIdxs[pt]])
            row+=velocityLocalBlockSize;
          else
          {
            typename FluidStateType::VelocityDiscreteFunctionType::RangeType temp;
            const auto xGlobal(entity.geometry().global(lagrangePointSet.point(pt)));
            #if INTERPOLATION_TYPE == 0
            oldFluidState.velocity().evaluate(xGlobal,temp);
            #else
            averageResearchDepth+=searchForEntity(oldFluidState.bulkGridPart(),std::move(oldEntity),xGlobal);
            auto localOldVelocity(oldFluidState.velocity().localFunction(oldEntity));
            localOldVelocity.evaluate(oldEntity.geometry().local(xGlobal),temp);
            #endif
            for(auto l=0;l!=velocityLocalBlockSize;++l,++row)
              localVelocity[row]=temp[l];
            dofAlreadyInterpolated[globalIdxs[pt]]=true;
          }
        }
      }
      #if INTERPOLATION_TYPE != 0
      std::cout<<"Average research depth : "<<static_cast<double>(averageResearchDepth)/static_cast<double>(dofAlreadyInterpolated.size())
        <<"."<<std::endl;
      #endif
    }
    timerInterpolation.stop();
    timerStep.stop();

    if(!(femScheme.problem().isDensityNull()))
      std::cout<<"Velocity interpolation time: "<<timerInterpolation.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Full timestep time: "<<timerStep.elapsed()<<" seconds."<<std::endl;
  }

  // print errors
  if(computeErrors)
  {
    errors[1]=std::pow(errors[1]*timeProvider.deltaT(),0.5);
    errors[2]=std::pow(errors[2]*timeProvider.deltaT(),0.5);
    errors[4]=std::pow(errors[4]*timeProvider.deltaT(),0.5);
    std::cout<<std::endl<<"Errors:"<<std::endl;
    std::cout.precision(5);
    std::cout<<std::scientific;
    #if PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 3
    std::cout<<"||X-x||_{L^oo} = "<<errors[0]<<std::endl;
    #endif
    std::cout<<"||U-Iu||_{L^2} = "<<errors[1]<<std::endl;
    std::cout<<"||U-Iu||_{H^1} = "<<errors[2]<<std::endl;
    std::cout<<"||U-Iu||_{L^oo} = "<<errors[3]<<std::endl;
    #if PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 3
    std::cout<<"||P-p||_{L^2} = "<<errors[4]<<std::endl;
    #endif
    std::cout<<"||P-Ip||_{L^oo} = "<<errors[5]<<std::endl;
    std::cout.unsetf(std::ios_base::floatfield);
  }
}

}
}

#endif // DUNE_FEM_COMPUTE_HH
