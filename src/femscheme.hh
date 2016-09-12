#ifndef DUEN_FEM_FEMSCHEME_HH
#define DUEN_FEM_FEMSCHEME_HH

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <dune/common/timer.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/istl/solvers.hh>
#include <dune/fem/solver/umfpacksolver.hh>
#include <dune/fem/solver/spqrsolver.hh>
#include <dune/fem/solver/ldlsolver.hh>

#include "problems.hh"
#include "nulloperator.hh"
#include "massmatrix.hh"
#include "operatorwrapper.hh"
#include "operatorgluer.hh"
#include "coupledoperatorwrapper.hh"
#include "preconditioners.hh"

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
  typedef typename FluidStateType::PressureDiscreteFunctionType PressureDiscreteFunctionType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef typename FluidStateType::PressureAdditionalDiscreteFunctionType PressureAdditionalDiscreteFunctionType;
  #endif
  typedef typename FluidStateType::CurvatureDiscreteFunctionType CurvatureDiscreteFunctionType;
  typedef typename FluidStateType::DisplacementDiscreteFunctionType DisplacementDiscreteFunctionType;
  typedef typename FluidStateType::BulkDiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename FluidStateType::InterfaceDiscreteFunctionType InterfaceDiscreteFunctionType;

  // define bulk-interface mapper
  typedef typename FluidStateType::CoupledMeshManagerType::BulkInterfaceGridMapperType BulkInterfaceGridMapperType;

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
  typedef NavierStokesTest1Problem<FluidStateType> ProblemType;
  #else
  typedef NavierStokesTest2Problem<FluidStateType> ProblemType;
  #endif

  // define operators for bulk
  typedef BulkVelocityOperator<VelocityDiscreteFunctionType,ProblemType> VelocityOperatorType;
  typedef BulkVelocityPressureOperator<VelocityDiscreteFunctionType,PressureDiscreteFunctionType> VelocityPressureOperatorType;
  typedef BulkPressureVelocityOperator<PressureDiscreteFunctionType,VelocityDiscreteFunctionType> PressureVelocityOperatorType;
  typedef NullOperator<PressureDiscreteFunctionType,PressureDiscreteFunctionType> PressureOperatorType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef BulkVelocityPressureOperator<VelocityDiscreteFunctionType,PressureAdditionalDiscreteFunctionType>
    VelocityPressureAdditionalOperatorType;
  typedef BulkPressureVelocityOperator<PressureAdditionalDiscreteFunctionType,VelocityDiscreteFunctionType>
    PressureAdditionalVelocityOperatorType;
  #endif
  typedef MassMatrix<PressureDiscreteFunctionType,PressureDiscreteFunctionType> BulkMassMatrixOperatorType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef MassMatrix<PressureAdditionalDiscreteFunctionType,PressureAdditionalDiscreteFunctionType> BulkMassMatrixAdditionalOperatorType;
  #endif
  typedef CurvatureVelocityOperator<CurvatureDiscreteFunctionType,VelocityDiscreteFunctionType,BulkInterfaceGridMapperType>
    CurvatureVelocityOperatorType;

  // define operators for interface
  typedef InterfaceOperator<InterfaceDiscreteFunctionType> InterfaceOperatorType;
  typedef VelocityCurvatureOperator<VelocityDiscreteFunctionType,CurvatureDiscreteFunctionType,BulkInterfaceGridMapperType>
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
    interfaceOp.assemble(fluidstate_.meshManager().mapper(),timeProvider);

    //assemble rhs
    InterfaceDiscreteFunctionType rhs("interface RHS",fluidstate_.interfaceSpace());
    rhs.clear();
    assembleInterfaceRHS(rhs,interfaceOp);
    rhs*=-1.0;

    // solve
    UMFPACKOp<InterfaceDiscreteFunctionType,InterfaceOperatorType> interfaceInvOp(interfaceOp);
    interfaceInvOp(rhs,fluidstate_.interfaceSolution());
  }

  // compute solution
  template<typename TimeProviderType>
  void operator()(const TimeProviderType& timeProvider)
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
    PressureVelocityOperatorType pressureVelocityOp(fluidstate_.pressureSpace(),fluidstate_.velocitySpace());
    pressureVelocityOp.assemble();
    #if PRESSURE_SPACE_TYPE == 2
    PressureAdditionalVelocityOperatorType pressureAdditionalVelocityOp(fluidstate_.pressureAdditionalSpace(),fluidstate_.velocitySpace());
    pressureAdditionalVelocityOp.assemble();
    #endif
    VelocityPressureOperatorType velocityPressureOp(fluidstate_.velocitySpace(),fluidstate_.pressureSpace());
    #if !USE_SYMMETRIC_DIRICHLET
    velocityPressureOp.assembleTransposingOp(pressureVelocityOp);
    #endif
    #if PRESSURE_SPACE_TYPE == 2
    VelocityPressureAdditionalOperatorType velocityPressureAdditionalOp(fluidstate_.velocitySpace(),fluidstate_.pressureAdditionalSpace());
    #if !USE_SYMMETRIC_DIRICHLET
    velocityPressureAdditionalOp.assembleTransposingOp(pressureAdditionalVelocityOp);
    #endif
    #endif
    PressureOperatorType pressureOp(fluidstate_.pressureSpace(),fluidstate_.pressureSpace());
    //pressureOp.assemble(); // not needed since it is never used
    CurvatureVelocityOperatorType curvatureVelocityOp(fluidstate_.curvatureSpace(),fluidstate_.velocitySpace(),
                                                      fluidstate_.meshManager().mapper());
    curvatureVelocityOp.assemble();
    BulkMassMatrixOperatorType bulkMassMatrixOp(fluidstate_.pressureSpace());
    bulkMassMatrixOp.assemble();
    #if PRESSURE_SPACE_TYPE == 2
    BulkMassMatrixAdditionalOperatorType bulkMassMatrixAdditionalOp(fluidstate_.pressureAdditionalSpace());
    bulkMassMatrixAdditionalOp.assemble();
    #endif
    timerAssembleBulk.stop();

    // assemble interface operators
    timerAssembleInterface.start();
    InterfaceOperatorType interfaceOp(fluidstate_.interfaceSpace());
    interfaceOp.assemble(fluidstate_.meshManager().mapper(),timeProvider);
    VelocityCurvatureOperatorType velocityCurvatureOp(fluidstate_.velocitySpace(),fluidstate_.curvatureSpace(),
                                                      fluidstate_.meshManager().mapper());
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
    assembleVelocityRHS(velocityRHS,fluidstate_.velocity(),problem_,timeProvider);
    #if PROBLEM_NUMBER == 3 || PROBLEM_NUMBER == 4 || PROBLEM_NUMBER == 8 || PROBLEM_NUMBER == 9
    assemblePressureRHS(bulkRHS.template subDiscreteFunction<1>(),problem_.velocityBC(),timeProvider);
    #if PRESSURE_SPACE_TYPE == 2
    assemblePressureRHS(bulkRHS.template subDiscreteFunction<2>(),problem_.velocityBC(),timeProvider);
    #endif
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
      velocityRHS.axpy(-1.0*gamma,velocityCouplingRHS);
      timerAssembleBulk.stop();
    }

    // impose bulk bc
    timerAssembleBulk.start();
    #if PRESSURE_SPACE_TYPE == 2
    problem_.applyBC(timeProvider,bulkRHS,velocityOp,pressureVelocityOp,pressureAdditionalVelocityOp);
    #else
    problem_.applyBC(timeProvider,bulkRHS,velocityOp,pressureVelocityOp);
    #endif
    timerAssembleBulk.stop();

    // assemble bulk operators
    #if USE_SYMMETRIC_DIRICHLET
    timerAssembleBulk.start();
    velocityPressureOp.assembleTransposingOp(pressureVelocityOp);
    #if PRESSURE_SPACE_TYPE == 2
    velocityPressureAdditionalOp.assembleTransposingOp(pressureAdditionalVelocityOp);
    #endif
    timerAssembleBulk.stop();
    #endif

    // compute bulk solution
    const int verbosity(Parameter::getValue<int>("SolverVerbosity",0));
    const int maxIter(Parameter::getValue<int>("SolverMaxIter",1000));
    const double redEps(Parameter::getValue<double>("SolverReductionEpsilon",1.e-12));
    const int restart(Parameter::getValue<int>("SolverRestart",5));

    timerSolveBulk.start();
    typedef CoupledOperatorWrapper<VelocityOperatorType,CurvatureVelocityOperatorType,InterfaceOperatorType,
                                   InterfaceInverseOperatorType,VelocityCurvatureOperatorType> CoupledOperatorWrapperType;
    CoupledOperatorWrapperType coupledWrapperOp(velocityOp,curvatureVelocityOp,interfaceOp,interfaceInvOp,velocityCurvatureOp,gamma);

    #if PRESSURE_SPACE_TYPE != 2
    typedef OperatorWrapper<CoupledOperatorWrapperType,PressureVelocityOperatorType,
                            VelocityPressureOperatorType,PressureOperatorType> BulkOperatorWrapperType;
    BulkOperatorWrapperType bulkOp(coupledWrapperOp,pressureVelocityOp,velocityPressureOp,pressureOp);
    #if PRECONDITIONER_TYPE == 0
    typedef StokesPrecond<BulkDiscreteFunctionType,VelocityOperatorType,PressureVelocityOperatorType,BulkMassMatrixOperatorType>
      BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(velocityOp,pressureVelocityOp,bulkMassMatrixOp);
    #elif PRECONDITIONER_TYPE == 1
    pressureOp.assemble();
    typedef OperatorGluer<VelocityOperatorType,PressureVelocityOperatorType,VelocityPressureOperatorType,PressureOperatorType>
      OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressureVelocityOp,velocityPressureOp,pressureOp);
    opGluer.assemble();
    typedef DirectPrecond<OperatorGluerType,UMFPACKOp> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #elif PRECONDITIONER_TYPE == 2
    pressureOp.assemble();
    typedef OperatorGluer<VelocityOperatorType,PressureVelocityOperatorType,VelocityPressureOperatorType,PressureOperatorType>
      OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressureVelocityOp,velocityPressureOp,pressureOp);
    opGluer.assemble();
    typedef DirectPrecond<OperatorGluerType,SPQROp> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #else
    pressureOp.assemble();
    typedef OperatorGluer<VelocityOperatorType,PressureVelocityOperatorType,VelocityPressureOperatorType,PressureOperatorType>
      OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressureVelocityOp,velocityPressureOp,pressureOp);
    opGluer.assemble();
    typedef DirectPrecond<OperatorGluerType,LDLOp> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #endif
    #else
    typedef ExtendedOperatorWrapper<CoupledOperatorWrapperType,PressureVelocityOperatorType,VelocityPressureOperatorType,
      PressureAdditionalVelocityOperatorType,VelocityPressureAdditionalOperatorType> BulkOperatorWrapperType;
    BulkOperatorWrapperType bulkOp(coupledWrapperOp,pressureVelocityOp,velocityPressureOp,pressureAdditionalVelocityOp,
                                   velocityPressureAdditionalOp);
    #if PRECONDITIONER_TYPE == 0
    typedef ExtendedStokesPrecond<BulkDiscreteFunctionType,VelocityOperatorType,PressureVelocityOperatorType,BulkMassMatrixOperatorType,
                                  PressureAdditionalVelocityOperatorType,BulkMassMatrixAdditionalOperatorType> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(velocityOp,pressureVelocityOp,bulkMassMatrixOp,pressureAdditionalVelocityOp,
                                              bulkMassMatrixAdditionalOp);
    #elif PRECONDITIONER_TYPE == 1
    typedef ExtendedOperatorGluer<VelocityOperatorType,PressureVelocityOperatorType,VelocityPressureOperatorType,
                                  PressureAdditionalVelocityOperatorType,VelocityPressureAdditionalOperatorType> OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressureVelocityOp,velocityPressureOp,pressureAdditionalVelocityOp,velocityPressureAdditionalOp);
    opGluer.assemble();
    typedef DirectPrecond<OperatorGluerType,UMFPACKOp> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #elif PRECONDITIONER_TYPE == 2
    typedef ExtendedOperatorGluer<VelocityOperatorType,PressureVelocityOperatorType,VelocityPressureOperatorType,
                                  PressureAdditionalVelocityOperatorType,VelocityPressureAdditionalOperatorType> OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressureVelocityOp,velocityPressureOp,pressureAdditionalVelocityOp,velocityPressureAdditionalOp);
    opGluer.assemble();
    typedef DirectPrecond<OperatorGluerType,SPQROp> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #else
    typedef ExtendedOperatorGluer<VelocityOperatorType,PressureVelocityOperatorType,VelocityPressureOperatorType,
                                  PressureAdditionalVelocityOperatorType,VelocityPressureAdditionalOperatorType> OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressureVelocityOp,velocityPressureOp,pressureAdditionalVelocityOp,velocityPressureAdditionalOp);
    opGluer.assemble();
    typedef DirectPrecond<OperatorGluerType,LDLOp> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #endif
    #endif

    InverseOperatorResult returnInfo;
    Dune::RestartedGMResSolver<BulkDiscreteFunctionType> bulkInvOp(bulkOp,bulkPreconditioner,redEps,restart,maxIter,verbosity);
    bulkInvOp.apply(fluidstate_.bulkSolution(),bulkRHS,returnInfo);
    timerSolveBulk.stop();

    // add interface coupling
    timerAssembleInterface.start();
    velocityCurvatureOp(fluidstate_.velocity(),interfaceRHS.template subDiscreteFunction<0>());
    interfaceRHS*=-1.0;
    timerAssembleInterface.stop();

    // solve interface
    timerSolveInterface.start();
    interfaceInvOp.apply(interfaceRHS,fluidstate_.interfaceSolution());
    interfaceInvOp.finalize();
    timerSolveInterface.stop();

    // project pressure solution to the space of mean zero function
    timerSolveBulk.start();
    projectZeroMean(fluidstate_.pressure(),bulkMassMatrixOp);
    #if PRESSURE_SPACE_TYPE == 2
    projectZeroMean(fluidstate_.pressureAdditional(),bulkMassMatrixAdditionalOp);
    #endif
    timerSolveBulk.stop();

    // remove displacement if the problem has only 1 phase
    if(gamma==0.0)
      fluidstate_.displacement().clear();

    // print timers
    std::cout<<"Assemble bulk operators (";
    fluidstate_.printBulkInfo();
    std::cout<<") time: "<<timerAssembleBulk.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Assemble interface operator (";
    fluidstate_.printInterfaceInfo();
    std::cout<<") time: "<<timerAssembleInterface.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Solve bulk time: "<<timerSolveBulk.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Solve interface time: "<<timerSolveInterface.elapsed()<<" seconds."<<std::endl;
  }

  private:
  FluidStateType& fluidstate_;
  ProblemType problem_;

  // project df into the space of mean zero function
  template<typename DF,typename Op>
  void projectZeroMean(DF& df,const Op& op) const
  {
    DF temp("temp",df.space());
    DF ones("ones",df.space());
    std::fill(ones.dbegin(),ones.dend(),1.0);
    op(df,temp);
    auto coeff(ones.scalarProductDofs(temp));
    op(ones,temp);
    coeff*=-1.0/ones.scalarProductDofs(temp);
    df.axpy(coeff,ones);
  }
};

}
}

#endif // DUEN_FEM_FEMSCHEME_HH
