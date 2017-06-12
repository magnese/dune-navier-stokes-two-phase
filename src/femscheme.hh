#ifndef DUEN_FEM_FEMSCHEME_HH
#define DUEN_FEM_FEMSCHEME_HH

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/timer.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/integrator.hh>
#include <dune/fem/solver/umfpacksolver.hh>
#include <dune/fem/solver/spqrsolver.hh>
#include <dune/fem/space/common/interpolate.hh>

#include "problems.hh"
#include "nulloperator.hh"
#include "massmatrix.hh"
#include "operatorwrapper.hh"
#include "operatorgluer.hh"
#include "coupledoperatorwrapper.hh"
#include "preconditioners.hh"
#include "restartedgmrescaller.hh"

#include "bulkvelocityoperator.hh"
#include "bulkvelocitypressureoperator.hh"
#include "bulkpressurevelocityoperator.hh"
#include "assemblevelocityrhs.hh"
#include "assemblepressurerhs.hh"

#include "interfaceoperator.hh"
#include "assembleinterfacerhs.hh"

#include "curvaturevelocityoperator.hh"
#include "velocitycurvatureoperator.hh"

namespace Dune
{
namespace Fem
{

template<typename FluidStateImp>
class FemScheme
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef FemScheme<FluidStateType> ThisType;

  // define discrete functions
  typedef typename FluidStateType::VelocityDiscreteFunctionType VelocityDiscreteFunctionType;
  typedef typename FluidStateType::Pressure0DiscreteFunctionType Pressure0DiscreteFunctionType;
  #if USE_EXTENDED_PRESSURE_SPACE
  typedef typename FluidStateType::Pressure1DiscreteFunctionType Pressure1DiscreteFunctionType;
  #endif
  typedef typename FluidStateType::CurvatureDiscreteFunctionType CurvatureDiscreteFunctionType;
  typedef typename FluidStateType::DisplacementDiscreteFunctionType DisplacementDiscreteFunctionType;
  typedef typename FluidStateType::BulkDiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename FluidStateType::InterfaceDiscreteFunctionType InterfaceDiscreteFunctionType;

  // define coupled mesh manager
  typedef typename FluidStateType::CoupledMeshManagerType CoupledMeshManagerType;

  // define problem
  #if PROBLEM_NUMBER == 0
  typedef StokesTest1Problem<FluidStateType> ProblemType;
  #elif PROBLEM_NUMBER == 1
  typedef StokesTest2Problem<FluidStateType> ProblemType;
  #elif PROBLEM_NUMBER == 2
  typedef StationaryBubbleProblem<FluidStateType> ProblemType;
  #elif PROBLEM_NUMBER == 3
  typedef ExpandingBubbleProblem<FluidStateType> ProblemType;
  #elif PROBLEM_NUMBER == 4
  typedef ShearFlowProblem<FluidStateType> ProblemType;
  #elif PROBLEM_NUMBER == 5
  typedef StationaryNavierStokesProblem<FluidStateType> ProblemType;
  #elif PROBLEM_NUMBER == 6
  typedef NavierStokes2DProblem<FluidStateType> ProblemType;
  #elif PROBLEM_NUMBER == 7
  typedef RisingBubbleProblem<FluidStateType> ProblemType;
  #elif PROBLEM_NUMBER == 8
  typedef NavierStokesExpandingBubble1Problem<FluidStateType> ProblemType;
  #elif PROBLEM_NUMBER == 9
  typedef NavierStokesExpandingBubble2Problem<FluidStateType> ProblemType;
  #elif PROBLEM_NUMBER == 10
  typedef NavierStokesExpandingBubble3Problem<FluidStateType> ProblemType;
  #endif

  // define operators for bulk
  typedef BulkVelocityOperator<VelocityDiscreteFunctionType,ProblemType> VelocityOperatorType;
  typedef BulkVelocityPressureOperator<VelocityDiscreteFunctionType,Pressure0DiscreteFunctionType> VelocityPressure0OperatorType;
  typedef BulkPressureVelocityOperator<Pressure0DiscreteFunctionType,VelocityDiscreteFunctionType> Pressure0VelocityOperatorType;
  typedef NullOperator<Pressure0DiscreteFunctionType,Pressure0DiscreteFunctionType> Pressure0OperatorType;
  #if USE_EXTENDED_PRESSURE_SPACE
  typedef BulkVelocityPressureOperator<VelocityDiscreteFunctionType,Pressure1DiscreteFunctionType> VelocityPressure1OperatorType;
  typedef BulkPressureVelocityOperator<Pressure1DiscreteFunctionType,VelocityDiscreteFunctionType> Pressure1VelocityOperatorType;
  #endif
  typedef MassMatrix<Pressure0DiscreteFunctionType,Pressure0DiscreteFunctionType> BulkMassMatrixOperatorType;
  #if USE_EXTENDED_PRESSURE_SPACE
  typedef MassMatrix<Pressure1DiscreteFunctionType,Pressure1DiscreteFunctionType> BulkMassMatrixAdditionalOperatorType;
  #endif
  typedef CurvatureVelocityOperator<CurvatureDiscreteFunctionType,VelocityDiscreteFunctionType,CoupledMeshManagerType>
    CurvatureVelocityOperatorType;

  // define operators for interface
  typedef InterfaceOperator<InterfaceDiscreteFunctionType> InterfaceOperatorType;
  typedef VelocityCurvatureOperator<VelocityDiscreteFunctionType,CurvatureDiscreteFunctionType,CoupledMeshManagerType>
    VelocityCurvatureOperatorType;

  // constructor
  explicit FemScheme(FluidStateType& fluidState):
    fluidstate_(fluidState),problem_(fluidState)
  {}

  FemScheme(const ThisType& )=delete;

  // get fluid state
  const FluidStateType& fluidState() const
  {
    return fluidstate_;
  }
  FluidStateType& fluidState()
  {
    return fluidstate_;
  }

  // get problem
  const ProblemType& problem() const
  {
    return problem_;
  }
  ProblemType& problem()
  {
    return problem_;
  }

  // compute interface initial curvature
  template<typename TimeProviderType>
  void computeInterface(const TimeProviderType& timeProvider)
  {
    // rebuild all quantities if the mesh is changed
    fluidstate_.update();

    // assemble operator
    InterfaceOperatorType interfaceOp(fluidstate_.interfaceSpace());
    interfaceOp.assemble(fluidstate_.meshManager(),timeProvider);

    // assemble rhs
    InterfaceDiscreteFunctionType rhs("interface RHS",fluidstate_.interfaceSpace());
    rhs.clear();
    assembleInterfaceRHS(rhs,interfaceOp);

    // solve
    UMFPACKOp<InterfaceDiscreteFunctionType,InterfaceOperatorType> interfaceInvOp(interfaceOp);
    interfaceInvOp(rhs,fluidstate_.interfaceSolution());
  }

  // compute solution
  template<typename TimeProviderType,typename MeshSmoothingType>
  void operator()(const TimeProviderType& timeProvider,MeshSmoothingType& meshSmoothing)
  {
    // rebuild all quantities if the mesh is changed
    fluidstate_.update();

    // create timers
    Timer timerAssembleBulk(false);
    Timer timerAssembleInterface(false);
    Timer timerSolveBulk(false);
    Timer timerSolveInterface(false);

    // assemble bulk operators
    timerAssembleBulk.start();
    VelocityOperatorType velocityOp(fluidstate_.velocitySpace(),problem_,fluidstate_.velocity());
    velocityOp.assemble(timeProvider);
    Pressure0VelocityOperatorType pressure0VelocityOp(fluidstate_.pressure0Space(),fluidstate_.velocitySpace());
    pressure0VelocityOp.assemble();
    #if USE_EXTENDED_PRESSURE_SPACE
    Pressure1VelocityOperatorType pressure1VelocityOp(fluidstate_.pressure1Space(),fluidstate_.velocitySpace());
    pressure1VelocityOp.assemble();
    #endif
    VelocityPressure0OperatorType velocityPressure0Op(fluidstate_.velocitySpace(),fluidstate_.pressure0Space());
    #if !USE_SYMMETRIC_DIRICHLET
    velocityPressure0Op.assembleTransposingOp(pressure0VelocityOp);
    #endif
    #if USE_EXTENDED_PRESSURE_SPACE
    VelocityPressure1OperatorType velocityPressure1Op(fluidstate_.velocitySpace(),fluidstate_.pressure1Space());
    #if !USE_SYMMETRIC_DIRICHLET
    velocityPressure1Op.assembleTransposingOp(pressure1VelocityOp);
    #endif
    #endif
    Pressure0OperatorType pressure0Op(fluidstate_.pressure0Space(),fluidstate_.pressure0Space());
    CurvatureVelocityOperatorType curvatureVelocityOp(fluidstate_.curvatureSpace(),fluidstate_.velocitySpace(),fluidstate_.meshManager());
    curvatureVelocityOp.assemble();
    BulkMassMatrixOperatorType bulkMassMatrixOp(fluidstate_.pressure0Space());
    bulkMassMatrixOp.assemble();
    #if USE_EXTENDED_PRESSURE_SPACE
    BulkMassMatrixAdditionalOperatorType bulkMassMatrixAdditionalOp(fluidstate_.pressure1Space());
    bulkMassMatrixAdditionalOp.assemble();
    #endif
    timerAssembleBulk.stop();

    // assemble interface operators
    timerAssembleInterface.start();
    InterfaceOperatorType interfaceOp(fluidstate_.interfaceSpace());
    interfaceOp.assemble(fluidstate_.meshManager(),timeProvider);
    VelocityCurvatureOperatorType velocityCurvatureOp(fluidstate_.velocitySpace(),fluidstate_.curvatureSpace(),fluidstate_.meshManager());
    velocityCurvatureOp.assemble();
    timerAssembleInterface.stop();

    // assemble interface inverse operator
    timerSolveInterface.start();
    typedef UMFPACKOp<InterfaceDiscreteFunctionType,InterfaceOperatorType> InterfaceInverseOperatorType;
    InterfaceInverseOperatorType interfaceInvOp(interfaceOp);
    interfaceInvOp.prepare();
    timerSolveInterface.stop();

    // assemble bulk RHS
    timerAssembleBulk.start();
    BulkDiscreteFunctionType bulkRHS("bulk RHS",fluidstate_.bulkSpace());
    bulkRHS.clear();
    auto& velocityRHS(bulkRHS.template subDiscreteFunction<0>());
    assembleVelocityRHS(velocityRHS,fluidstate_,problem_,timeProvider);
    #if PROBLEM_NUMBER == 3 || PROBLEM_NUMBER == 4 || PROBLEM_NUMBER == 8 || PROBLEM_NUMBER == 9 || PROBLEM_NUMBER == 10
    Hybrid::forEach(std::make_index_sequence<BulkDiscreteFunctionType::Sequence::size()-1>{},
      [&](auto i){assemblePressureRHS(bulkRHS.template subDiscreteFunction<i+1>(),problem_.velocityBC(),timeProvider);});
    #endif
    timerAssembleBulk.stop();

    // assemble interface RHS
    timerAssembleInterface.start();
    InterfaceDiscreteFunctionType interfaceRHS("interface RHS",fluidstate_.interfaceSpace());
    interfaceRHS.clear();
    assembleInterfaceRHS(interfaceRHS,interfaceOp);
    timerAssembleInterface.stop();

    // add bulk coupling
    const auto gamma(problem_.gamma());
    if(gamma!=0.0)
    {
      timerAssembleBulk.start();
      InterfaceDiscreteFunctionType interfaceTempFunction("interface temporary function",fluidstate_.interfaceSpace());
      interfaceInvOp.apply(interfaceRHS,interfaceTempFunction);
      VelocityDiscreteFunctionType velocityCouplingRHS("velocity coupling RHS",fluidstate_.velocitySpace());
      curvatureVelocityOp(interfaceTempFunction.template subDiscreteFunction<0>(),velocityCouplingRHS);
      velocityRHS.axpy(gamma,velocityCouplingRHS);
      timerAssembleBulk.stop();
    }

    // impose bulk bc
    timerAssembleBulk.start();
    #if !USE_EXTENDED_PRESSURE_SPACE
    problem_.applyBC(timeProvider,bulkRHS,velocityOp,pressure0VelocityOp);
    #else
    problem_.applyBC(timeProvider,bulkRHS,velocityOp,pressure0VelocityOp,pressure1VelocityOp);
    #endif
    timerAssembleBulk.stop();

    // assemble bulk operators
    #if USE_SYMMETRIC_DIRICHLET
    timerAssembleBulk.start();
    velocityPressure0Op.assembleTransposingOp(pressure0VelocityOp);
    #if USE_EXTENDED_PRESSURE_SPACE
    velocityPressure1Op.assembleTransposingOp(pressure1VelocityOp);
    #endif
    timerAssembleBulk.stop();
    #endif

    // create operators wrappers and preconditioners
    timerSolveBulk.start();
    typedef CoupledOperatorWrapper<VelocityOperatorType,CurvatureVelocityOperatorType,InterfaceOperatorType,
                                   InterfaceInverseOperatorType,VelocityCurvatureOperatorType> CoupledOperatorWrapperType;
    CoupledOperatorWrapperType coupledWrapperOp(velocityOp,curvatureVelocityOp,interfaceOp,interfaceInvOp,velocityCurvatureOp,gamma);

    #if !USE_EXTENDED_PRESSURE_SPACE
    typedef OperatorWrapper<CoupledOperatorWrapperType,Pressure0VelocityOperatorType,
                            VelocityPressure0OperatorType,Pressure0OperatorType> BulkOperatorWrapperType;
    BulkOperatorWrapperType bulkOp(coupledWrapperOp,pressure0VelocityOp,velocityPressure0Op,pressure0Op);
    #if PRECONDITIONER_TYPE == 0
    typedef StokesPrecond<BulkDiscreteFunctionType,VelocityOperatorType,Pressure0VelocityOperatorType,BulkMassMatrixOperatorType>
      BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(velocityOp,pressure0VelocityOp,bulkMassMatrixOp);
    #else
    typedef OperatorGluer<VelocityOperatorType,Pressure0VelocityOperatorType,VelocityPressure0OperatorType,Pressure0OperatorType>
      OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressure0VelocityOp,velocityPressure0Op,pressure0Op);
    opGluer.assemble();
    #if PRECONDITIONER_TYPE == 1
    typedef DirectPrecond<OperatorGluerType,UMFPACKOp> BulkPreconditionerType;
    #elif PRECONDITIONER_TYPE == 2
    typedef DirectPrecond<OperatorGluerType,SPQROp> BulkPreconditionerType;
    #endif
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #endif
    #else
    typedef ExtendedOperatorWrapper<CoupledOperatorWrapperType,Pressure0VelocityOperatorType,VelocityPressure0OperatorType,
      Pressure1VelocityOperatorType,VelocityPressure1OperatorType> BulkOperatorWrapperType;
    BulkOperatorWrapperType bulkOp(coupledWrapperOp,pressure0VelocityOp,velocityPressure0Op,pressure1VelocityOp,
                                   velocityPressure1Op);
    #if PRECONDITIONER_TYPE == 0
    typedef ExtendedStokesPrecond<BulkDiscreteFunctionType,VelocityOperatorType,Pressure0VelocityOperatorType,BulkMassMatrixOperatorType,
                                  Pressure1VelocityOperatorType,BulkMassMatrixAdditionalOperatorType> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(velocityOp,pressure0VelocityOp,bulkMassMatrixOp,pressure1VelocityOp,
                                              bulkMassMatrixAdditionalOp);
    #else
    typedef ExtendedOperatorGluer<VelocityOperatorType,Pressure0VelocityOperatorType,VelocityPressure0OperatorType,
                                  Pressure1VelocityOperatorType,VelocityPressure1OperatorType> OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressure0VelocityOp,velocityPressure0Op,pressure1VelocityOp,velocityPressure1Op);
    opGluer.assemble();
    #if PRECONDITIONER_TYPE == 1
    typedef DirectPrecond<OperatorGluerType,UMFPACKOp> BulkPreconditionerType;
    #elif PRECONDITIONER_TYPE == 2
    typedef DirectPrecond<OperatorGluerType,SPQROp> BulkPreconditionerType;
    #endif
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #endif
    #endif
    timerSolveBulk.stop();

    //solve bulk
    timerSolveBulk.start();
    RestartedGMResCaller<BulkDiscreteFunctionType> bulkInvOp(bulkOp,bulkPreconditioner);
    bulkInvOp(fluidstate_.bulkSolution(),bulkRHS);
    timerSolveBulk.stop();

    // fixed point iteration
    const int nonLinearSolverType(Parameter::getValue<int>("NonLinearSolverType",0));
    const bool useALE(nonLinearSolverType>1);
    if((!problem_.isDensityNull())&&(nonLinearSolverType>0))
    {
      #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
      DUNE_THROW(InvalidStateException,"ERROR: cannot use antisymmetric convective term with fixed point iteration!");
      #endif
      // read non-linear solver parameters
      const int nonLinearSolverVerbosity(Parameter::getValue<int>("NonLinearSolverVerbosity",0));
      const int nonLinearSolverMaxIterations(Parameter::getValue<int>("NonLinearSolverMaxIterations",1000));
      const double nonLinearSolverTolerance(Parameter::getValue<double>("NonLinearSolverTolerance",1.e-8));
      if(nonLinearSolverVerbosity>1)
        std::cout<<"Entering in the non-linear solver"<<(useALE?" (with ALE)":"")<<"\n";
      // force smoothing in ALE
      if(useALE)
        meshSmoothing.enable();
      // create a copy of the velocity and of the bulk displacement
      auto oldVelocity(fluidstate_.velocity());
      auto oldBulkDisplacement(fluidstate_.bulkDisplacement());
      oldBulkDisplacement.clear();
      // compute bulk solution
      bool doIteration(true);
      int iterationNumber(0);
      while(doIteration&&(iterationNumber<nonLinearSolverMaxIterations))
      {
        doIteration=false;
        ++iterationNumber;
        if(useALE)
        {
          // add interface coupling
          timerAssembleInterface.start();
          velocityCurvatureOp(fluidstate_.velocity(),interfaceRHS.template subDiscreteFunction<0>());
          interfaceRHS.template subDiscreteFunction<0>()*=-1.0;
          timerAssembleInterface.stop();
          // solve interface
          timerSolveInterface.start();
          interfaceInvOp.apply(interfaceRHS,fluidstate_.interfaceSolution());
          timerSolveInterface.stop();
          // compute smoothing, assemble time derivative on the moved grid and restore the grid
          fluidstate_.interfaceGrid().coordFunction()+=fluidstate_.displacement();
          meshSmoothing.computeBulkDisplacement();
          fluidstate_.bulkGrid().coordFunction()+=fluidstate_.bulkDisplacement();
          timerAssembleBulk.start();
          velocityOp.allocateAndAssembleTimeDerivative(timeProvider);
          timerAssembleBulk.stop();
          fluidstate_.interfaceGrid().coordFunction()-=fluidstate_.displacement();
          fluidstate_.bulkGrid().coordFunction()-=fluidstate_.bulkDisplacement();
          // add ALE contribution
          timerAssembleBulk.start();
          auto bulkDisplacementVelocity(fluidstate_.velocity());
          interpolate(fluidstate_.bulkDisplacement(),bulkDisplacementVelocity);
          bulkDisplacementVelocity/=timeProvider.deltaT();
          fluidstate_.velocity()-=bulkDisplacementVelocity;
          // re-assemble velocity operator
          velocityOp.assembleRemainingTerms(timeProvider);
          timerAssembleBulk.stop();
        }
        else
        {
          // re-assemble velocity operator
          timerAssembleBulk.start();
          velocityOp.assemble(timeProvider);
          timerAssembleBulk.stop();
        }
        // re-impose bulk bc
        timerAssembleBulk.start();
        #if !USE_EXTENDED_PRESSURE_SPACE
        problem_.applyBC(timeProvider,bulkRHS,velocityOp,pressure0VelocityOp);
        #else
        problem_.applyBC(timeProvider,bulkRHS,velocityOp,pressure0VelocityOp,pressure1VelocityOp);
        #endif
        timerAssembleBulk.stop();
        // solve bulk
        timerSolveBulk.start();
        bulkInvOp(fluidstate_.bulkSolution(),bulkRHS);
        timerSolveBulk.stop();
        // check if another iteration is needed
        double nonLinearSolverResidual(0);
        double ALEResidual(0);
        timerSolveBulk.start();
        // compute ||U-U_old||_Loo
        auto oldVelocityIt(oldVelocity.dbegin());
        for(const auto& dof:dofs(fluidstate_.velocity()))
          nonLinearSolverResidual=std::max(nonLinearSolverResidual,std::abs(dof-*(oldVelocityIt++)));
        // set old velocity to new velocity
        oldVelocity.assign(fluidstate_.velocity());
        if(nonLinearSolverResidual>nonLinearSolverTolerance)
          doIteration=true;
        // compute ||X_bulk-X_bulk_old||_Loo
        if(useALE)
        {
          auto oldBulkDisplacementIt(oldBulkDisplacement.dbegin());
          for(const auto& dof:dofs(fluidstate_.bulkDisplacement()))
            ALEResidual=std::max(ALEResidual,std::abs(dof-*(oldBulkDisplacementIt++)));
          // set old bulk displacement to new bulk displacement
          oldBulkDisplacement.assign(fluidstate_.bulkDisplacement());
          if(ALEResidual>nonLinearSolverTolerance)
            doIteration=true;
        }
        timerSolveBulk.stop();
        // output residuals
        if(doIteration&&(nonLinearSolverVerbosity>1))
        {
          std::cout<<"Iteration "<<iterationNumber<<" --> fixed-point residual = "<<nonLinearSolverResidual;
          if(useALE)
            std::cout<<"; ALE residual = "<<ALEResidual;
          std::cout<<"\n";
        }
        if((!doIteration)&&(nonLinearSolverVerbosity>0))
        {
          std::cout<<"Scheme converged to the solution with fixed-point residual "<<nonLinearSolverResidual;
          if(useALE)
            std::cout<<" and ALE residual "<<ALEResidual;
          std::cout<<" in "<<iterationNumber<<" iterations\n";
        }
      }
    }

    if(!useALE)
    {
      // add interface coupling
      timerAssembleInterface.start();
      velocityCurvatureOp(fluidstate_.velocity(),interfaceRHS.template subDiscreteFunction<0>());
      interfaceRHS.template subDiscreteFunction<0>()*=-1.0;
      timerAssembleInterface.stop();
      // solve interface
      timerSolveInterface.start();
      interfaceInvOp.apply(interfaceRHS,fluidstate_.interfaceSolution());
      timerSolveInterface.stop();
    }

    // finalize interface inverse operator
    timerSolveInterface.start();
    interfaceInvOp.finalize();
    timerSolveInterface.stop();

    // set pressure to correct values (if needed)
    #if PRESSURE_SPACE_TYPE == 2
    timerSolveBulk.start();
    for(const auto& entity:entities(fluidstate_.pressure0()))
    {
      auto localPressure0(fluidstate_.pressure0().localFunction(entity));
      auto localPressure1(fluidstate_.pressure1().localFunction(entity));
      auto localPressure(fluidstate_.pressure().localFunction(entity));
      for(auto i=decltype(localPressure.size()){0};i!=localPressure.size();++i)
        localPressure[i]=localPressure0[i]+localPressure1[0];
    }
    timerSolveBulk.stop();
    #elif PRESSURE_SPACE_TYPE == 3
    timerSolveBulk.start();
    std::vector<typename FluidStateType::PressureDiscreteSpaceType::RangeFieldType> localDOFs;
    localDOFs.reserve(fluidstate_.pressureSpace().blockMapper().maxNumDofs()*FluidStateType::PressureDiscreteSpaceType::localBlockSize);
    for(const auto& entity:fluidstate_.pressureSpace())
    {
      const auto interpolation(fluidstate_.pressureSpace().interpolation(entity));
      localDOFs.resize(fluidstate_.pressureSpace().basisFunctionSet(entity).size());
      if(fluidstate_.bulkInnerGridPart().contains(entity))
      {
        const auto localPressure(fluidstate_.pressure0().localFunction(entity));
        interpolation(localPressure,localDOFs);
      }
      else
      {
        const auto localPressure(fluidstate_.pressure1().localFunction(entity));
        interpolation(localPressure,localDOFs);
      }
      fluidstate_.pressure().setLocalDofs(entity,localDOFs);
    }
    timerSolveBulk.stop();
    #endif

    // project pressure solution to the space of mean zero function
    timerSolveBulk.start();
    Integrator<CachingQuadrature<typename FluidStateType::BulkGridPartType,0>> integrator(2*fluidstate_.pressureSpace().order()+1);
    typename FluidStateType::PressureDiscreteFunctionType::RangeType pressureIntegral(0);
    for(const auto& entity:fluidstate_.pressureSpace())
    {
      auto localPressure(fluidstate_.pressure().localFunction(entity));
      integrator.integrateAdd(entity,localPressure,pressureIntegral);
    }
    pressureIntegral/=fluidstate_.meshManager().bulkVolume();
    for(auto& dof:dofs(fluidstate_.pressure()))
      dof-=pressureIntegral;
    timerSolveBulk.stop();

    // remove displacement if the problem has only 1 phase
    if(gamma==0.0)
      fluidstate_.displacement().clear();

    // print timers
    std::cout<<"Assemble bulk operators (";
    fluidstate_.printBulkInfo();
    std::cout<<") time: "<<timerAssembleBulk.elapsed()<<" seconds.\n";
    std::cout<<"Assemble interface operator (";
    fluidstate_.printInterfaceInfo();
    std::cout<<") time: "<<timerAssembleInterface.elapsed()<<" seconds.\n";
    std::cout<<"Solve bulk time: "<<timerSolveBulk.elapsed()<<" seconds.\n";
    std::cout<<"Solve interface time: "<<timerSolveInterface.elapsed()<<" seconds.\n";
  }

  private:
  FluidStateType& fluidstate_;
  ProblemType problem_;
};

}
}

#endif // DUEN_FEM_FEMSCHEME_HH
