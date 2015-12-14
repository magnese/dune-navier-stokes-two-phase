#ifndef DUEN_FEM_FEMSCHEME_HH
#define DUEN_FEM_FEMSCHEME_HH

// C++ includes
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// dune includes
#include <dune/common/timer.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/istl/solvers.hh>
#include <dune/fem/solver/umfpacksolver.hh>
#include <dune/fem/solver/spqrsolver.hh>

// local includes
#include "problems.hh"
#include "nulloperator.hh"
#include "nullrhs.hh"
#include "massmatrix.hh"
#include "operatorwrapper.hh"
#include "operatorgluer.hh"
#include "coupledoperatorwrapper.hh"
#include "preconditioners.hh"

#include "bulkvelocityoperator.hh"
#include "bulkvelocitypressureoperator.hh"
#include "bulkpressurevelocityoperator.hh"
#include "bulkvelocityrhs.hh"
#include "bulkpressurerhs.hh"

#include "interfaceoperator.hh"
#include "interfacerhs.hh"

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

  // define grid parts
  typedef typename FluidStateType::BulkGridType BulkGridType;
  typedef typename FluidStateType::BulkGridPartType BulkGridPartType;
  typedef typename FluidStateType::InterfaceGridType InterfaceGridType;
  typedef typename FluidStateType::InterfaceGridPartType InterfaceGridPartType;

  // define spaces
  typedef typename FluidStateType::VelocityDiscreteSpaceType VelocityDiscreteSpaceType;
  typedef typename FluidStateType::PressureDiscreteSpaceType PressureDiscreteSpaceType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef typename FluidStateType::PressureAdditionalDiscreteSpaceType PressureAdditionalDiscreteSpaceType;
  #endif
  typedef typename FluidStateType::PressureDumpDiscreteSpaceType PressureDumpDiscreteSpaceType;
  typedef typename FluidStateType::CurvatureDiscreteSpaceType CurvatureDiscreteSpaceType;
  typedef typename FluidStateType::DisplacementDiscreteSpaceType DisplacementDiscreteSpaceType;
  typedef typename FluidStateType::BulkDiscreteSpaceType BulkDiscreteSpaceType;
  typedef typename FluidStateType::InterfaceDiscreteSpaceType InterfaceDiscreteSpaceType;

  // define discrete functions
  typedef typename FluidStateType::VelocityDiscreteFunctionType VelocityDiscreteFunctionType;
  typedef typename FluidStateType::PressureDiscreteFunctionType PressureDiscreteFunctionType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef typename FluidStateType::PressureAdditionalDiscreteFunctionType PressureAdditionalDiscreteFunctionType;
  #endif
  typedef typename FluidStateType::PressureDumpDiscreteFunctionType PressureDumpDiscreteFunctionType;
  typedef typename FluidStateType::CurvatureDiscreteFunctionType CurvatureDiscreteFunctionType;
  typedef typename FluidStateType::DisplacementDiscreteFunctionType DisplacementDiscreteFunctionType;
  typedef typename FluidStateType::BulkDiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename FluidStateType::InterfaceDiscreteFunctionType InterfaceDiscreteFunctionType;

  // define coupled mesh manager and bulk-interface mapper
  typedef typename FluidStateType::CoupledMeshManagerType CoupledMeshManagerType;
  typedef typename CoupledMeshManagerType::BulkInterfaceGridMapperType BulkInterfaceGridMapperType;

  // define problem
  #if PROBLEM_NUMBER == 0
  typedef StokesTest1Problem<VelocityDiscreteSpaceType,PressureDumpDiscreteSpaceType,CoupledMeshManagerType> ProblemType;
  #elif PROBLEM_NUMBER == 1
  typedef StokesTest2Problem<VelocityDiscreteSpaceType,PressureDumpDiscreteSpaceType,CoupledMeshManagerType> ProblemType;
  #elif PROBLEM_NUMBER == 2
  typedef StationaryBubbleProblem<VelocityDiscreteSpaceType,PressureDumpDiscreteSpaceType,CoupledMeshManagerType> ProblemType;
  #elif PROBLEM_NUMBER == 3
  typedef ExpandingBubbleProblem<VelocityDiscreteSpaceType,PressureDumpDiscreteSpaceType,CoupledMeshManagerType> ProblemType;
  #elif PROBLEM_NUMBER == 4
  typedef ShearFlowProblem<VelocityDiscreteSpaceType,PressureDumpDiscreteSpaceType,CoupledMeshManagerType> ProblemType;
  #elif PROBLEM_NUMBER == 5
  typedef StationaryNavierStokesProblem<VelocityDiscreteSpaceType,PressureDumpDiscreteSpaceType,CoupledMeshManagerType> ProblemType;
  #elif PROBLEM_NUMBER == 6
  typedef NavierStokes2DProblem<VelocityDiscreteSpaceType,PressureDumpDiscreteSpaceType,CoupledMeshManagerType> ProblemType;
  #else
  typedef RisingBubble2DProblem<VelocityDiscreteSpaceType,PressureDumpDiscreteSpaceType,CoupledMeshManagerType> ProblemType;
  #endif

  // define underlying matrices structures for bulk
  typedef SparseRowLinearOperator<VelocityDiscreteFunctionType,VelocityDiscreteFunctionType> VelocityLinearOperatorType;
  typedef SparseRowLinearOperator<VelocityDiscreteFunctionType,PressureDiscreteFunctionType> VelocityPressureLinearOperatorType;
  typedef SparseRowLinearOperator<PressureDiscreteFunctionType,VelocityDiscreteFunctionType> PressureVelocityLinearOpearatorType;
  typedef SparseRowLinearOperator<PressureDiscreteFunctionType,PressureDiscreteFunctionType> PressureLinearOperatorType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef SparseRowLinearOperator<VelocityDiscreteFunctionType,PressureAdditionalDiscreteFunctionType>
    VelocityPressureAdditionalLinearOperatorType;
  typedef SparseRowLinearOperator<PressureAdditionalDiscreteFunctionType,VelocityDiscreteFunctionType>
    PressureAdditionalVelocityLinearOpearatorType;
  typedef SparseRowLinearOperator<PressureAdditionalDiscreteFunctionType,PressureAdditionalDiscreteFunctionType>
    PressureAdditionalLinearOperatorType;
  #endif
  typedef SparseRowLinearOperator<CurvatureDiscreteFunctionType,VelocityDiscreteFunctionType> CurvatureVelocityLinearOperatorType;

  // define underlying matrices structures for interface
  typedef SparseRowLinearOperator<VelocityDiscreteFunctionType,CurvatureDiscreteFunctionType> VelocityCurvatureLinearOperatorType;
  typedef SparseRowLinearOperator<InterfaceDiscreteFunctionType,InterfaceDiscreteFunctionType> InterfaceLinearOperatorType;

  // define operators for bulk
  typedef BulkVelocityOperator<VelocityLinearOperatorType,ProblemType> VelocityOperatorType;
  typedef BulkVelocityPressureOperator<VelocityPressureLinearOperatorType> VelocityPressureOperatorType;
  typedef BulkPressureVelocityOperator<PressureVelocityLinearOpearatorType> PressureVelocityOperatorType;
  typedef NullOperator<PressureLinearOperatorType> PressureOperatorType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef BulkVelocityPressureOperator<VelocityPressureAdditionalLinearOperatorType> VelocityPressureAdditionalOperatorType;
  typedef BulkPressureVelocityOperator<PressureAdditionalVelocityLinearOpearatorType> PressureAdditionalVelocityOperatorType;
  #endif
  typedef MassMatrix<PressureLinearOperatorType> BulkMassMatrixOperatorType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef MassMatrix<PressureAdditionalLinearOperatorType> BulkMassMatrixAdditionalOperatorType;
  #endif
  typedef CurvatureVelocityOperator<CurvatureVelocityLinearOperatorType,BulkInterfaceGridMapperType> CurvatureVelocityOperatorType;
  typedef BulkVelocityRHS<VelocityDiscreteFunctionType,ProblemType> VelocityRHSType;
  #if PROBLEM_NUMBER == 0 || PROBLEM_NUMBER == 1 || PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 5 || PROBLEM_NUMBER == 6 || PROBLEM_NUMBER == 7
  typedef NullRHS<PressureDiscreteFunctionType> PressureRHSType;
  #else
  typedef BulkPressureRHS<PressureDiscreteFunctionType> PressureRHSType;
  #endif

  // define operators for interface
  typedef InterfaceOperator<InterfaceLinearOperatorType> InterfaceOperatorType;
  typedef VelocityCurvatureOperator<VelocityCurvatureLinearOperatorType,BulkInterfaceGridMapperType> VelocityCurvatureOperatorType;
  typedef InterfaceRHS<InterfaceDiscreteFunctionType> InterfaceRHSType;

  // constructor
  explicit FemScheme(FluidStateType& fluidState):
    fluidstate_(fluidState),problem_(fluidstate_.meshManager())
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
    InterfaceRHSType interfaceRHS(fluidstate_.interfaceSpace());
    interfaceRHS.assemble(interfaceOp);
    interfaceRHS.rhs()*=-1.0;

    // solve
    typedef UMFPACKOp<InterfaceDiscreteFunctionType,InterfaceOperatorType> InterfaceInverseOperatorType;
    InterfaceInverseOperatorType interfaceInvOp(interfaceOp);
    interfaceInvOp(interfaceRHS.rhs(),fluidstate_.interfaceSolution());

    // set the fluid state for the interface with the correct quantities
    fluidstate_.finalizeInterfaceQuantities();
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

    // assemble bulk operators and impose BC
    timerAssembleBulk.start();
    VelocityOperatorType velocityOp(fluidstate_.velocitySpace(),problem_,fluidstate_.velocity());
    velocityOp.assemble(timeProvider);
    PressureVelocityOperatorType pressureVelocityOp(fluidstate_.pressureSpace(),fluidstate_.velocitySpace());
    pressureVelocityOp.assemble();
    #if PRESSURE_SPACE_TYPE == 2
    PressureAdditionalVelocityOperatorType pressureAdditionalVelocityOp(fluidstate_.pressureAdditionalSpace(),fluidstate_.velocitySpace());
    pressureAdditionalVelocityOp.assemble();
    problem_.applyBCToOperator(velocityOp,pressureVelocityOp,pressureAdditionalVelocityOp);
    #else
    problem_.applyBCToOperator(velocityOp,pressureVelocityOp);
    #endif
    VelocityPressureOperatorType velocityPressureOp(fluidstate_.velocitySpace(),fluidstate_.pressureSpace());
    velocityPressureOp.assemble();
    #if PRESSURE_SPACE_TYPE == 2
    VelocityPressureAdditionalOperatorType velocityPressureAdditionalOp(fluidstate_.velocitySpace(),fluidstate_.pressureAdditionalSpace());
    velocityPressureAdditionalOp.assemble();
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
    VelocityRHSType velocityRHS(fluidstate_.velocitySpace(),problem_,fluidstate_.velocity());
    velocityRHS.assemble(problem_.velocityRHS(),timeProvider);
    PressureRHSType pressureRHS(fluidstate_.pressureSpace());
    pressureRHS.assemble(problem_.velocityBC(),timeProvider);
    timerAssembleBulk.stop();

    // assemble interface RHS
    timerAssembleInterface.start();
    InterfaceRHSType interfaceRHS(fluidstate_.interfaceSpace());
    interfaceRHS.assemble(interfaceOp);
    timerAssembleInterface.stop();

    // add bulk coupling
    const auto gamma(problem_.gamma());
    if(gamma!=0.0)
    {
      timerAssembleBulk.start();
      InterfaceDiscreteFunctionType interfaceTempFunction("interface temporary function",fluidstate_.interfaceSpace());
      interfaceInvOp.apply(interfaceRHS.rhs(),interfaceTempFunction);

      CurvatureDiscreteFunctionType curvatureTempFunction("curvature temporary function",fluidstate_.curvatureSpace());
      std::copy_n(interfaceTempFunction.dbegin(),curvatureTempFunction.size(),curvatureTempFunction.dbegin());

      VelocityDiscreteFunctionType velocityCouplingRHS("velocity coupling RHS",fluidstate_.velocitySpace());
      curvatureVelocityOp(curvatureTempFunction,velocityCouplingRHS);
      velocityRHS.rhs().axpy(-1.0*gamma,velocityCouplingRHS);
      timerAssembleBulk.stop();
    }

    // impose bulk bc
    timerAssembleBulk.start();
    problem_.applyBCToRHS(velocityRHS.rhs(),timeProvider);
    timerAssembleBulk.stop();

    // create combined bulk RHS and copy RHS in the combined function
    timerSolveBulk.start();
    BulkDiscreteFunctionType bulkRHS("bulk RHS",fluidstate_.bulkSpace());
    auto bulkIt=std::copy(velocityRHS.rhs().dbegin(),velocityRHS.rhs().dend(),bulkRHS.dbegin());
    std::copy(pressureRHS.rhs().dbegin(),pressureRHS.rhs().dend(),bulkIt);
    timerSolveBulk.stop();

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
    typedef StokesPrecond<VelocityOperatorType,PressureVelocityOperatorType,BulkMassMatrixOperatorType> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(velocityOp,pressureVelocityOp,bulkMassMatrixOp);
    #elif PRECONDITIONER_TYPE == 1
    pressureOp.assemble();
    typedef OperatorGluer<VelocityOperatorType,PressureVelocityOperatorType,VelocityPressureOperatorType,PressureOperatorType>
      OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressureVelocityOp,velocityPressureOp,pressureOp);
    opGluer.assemble();
    opGluer.applyDoctoring();
    typedef DirectPrecond<OperatorGluerType,UMFPACKOp> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #else
    pressureOp.assemble();
    typedef OperatorGluer<VelocityOperatorType,PressureVelocityOperatorType,VelocityPressureOperatorType,PressureOperatorType>
      OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressureVelocityOp,velocityPressureOp,pressureOp);
    opGluer.assemble();
    typedef DirectPrecond<OperatorGluerType,SPQROp> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #endif
    #else
    typedef ExtendedOperatorWrapper<CoupledOperatorWrapperType,PressureVelocityOperatorType,VelocityPressureOperatorType,
                                    PressureAdditionalVelocityOperatorType,VelocityPressureAdditionalOperatorType> BulkOperatorWrapperType;
    BulkOperatorWrapperType bulkOp(coupledWrapperOp,pressureVelocityOp,velocityPressureOp,pressureAdditionalVelocityOp,
                                   velocityPressureAdditionalOp);
    #if PRECONDITIONER_TYPE == 0
    typedef ExtendedStokesPrecond<VelocityOperatorType,PressureVelocityOperatorType,BulkMassMatrixOperatorType,
                                  PressureAdditionalVelocityOperatorType,BulkMassMatrixAdditionalOperatorType> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(velocityOp,pressureVelocityOp,bulkMassMatrixOp,pressureAdditionalVelocityOp,
                                              bulkMassMatrixAdditionalOp);
    #elif PRECONDITIONER_TYPE == 1
    typedef ExtendedOperatorGluer<VelocityOperatorType,PressureVelocityOperatorType,VelocityPressureOperatorType,
                                  PressureAdditionalVelocityOperatorType,VelocityPressureAdditionalOperatorType> OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressureVelocityOp,velocityPressureOp,pressureAdditionalVelocityOp,velocityPressureAdditionalOp);
    opGluer.assemble();
    opGluer.applyDoctoring();
    typedef DirectPrecond<OperatorGluerType,UMFPACKOp> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #else
    typedef ExtendedOperatorGluer<VelocityOperatorType,PressureVelocityOperatorType,VelocityPressureOperatorType,
                                  PressureAdditionalVelocityOperatorType,VelocityPressureAdditionalOperatorType> OperatorGluerType;
    OperatorGluerType opGluer(velocityOp,pressureVelocityOp,velocityPressureOp,pressureAdditionalVelocityOp,velocityPressureAdditionalOp);
    opGluer.assemble();
    typedef DirectPrecond<OperatorGluerType,SPQROp> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(opGluer);
    #endif
    #endif

    typedef Dune::RestartedGMResSolver<typename BulkDiscreteFunctionType::DofStorageType> BulkLinearInverseOperatorType;
    InverseOperatorResult returnInfo;
    BulkLinearInverseOperatorType bulkInvOp(bulkOp,bulkPreconditioner,redEps,restart,maxIter,verbosity);
    bulkInvOp.apply(fluidstate_.bulkSolution().blockVector(),bulkRHS.blockVector(),returnInfo);
    timerSolveBulk.stop();

    // set the fluid state for the bulk with the correct quantities
    fluidstate_.finalizeBulkQuantities();
    timerSolveBulk.stop();

    // add interface coupling
    timerAssembleInterface.start();
    CurvatureDiscreteFunctionType curvatureCouplingRHS("curvature coupling RHS",fluidstate_.curvatureSpace());
    velocityCurvatureOp(fluidstate_.velocity(),curvatureCouplingRHS);
    std::copy(curvatureCouplingRHS.dbegin(),curvatureCouplingRHS.dend(),interfaceRHS.rhs().dbegin());
    interfaceRHS.rhs()*=-1.0;
    timerAssembleInterface.stop();

    // solve interface
    timerSolveInterface.start();
    interfaceInvOp.apply(interfaceRHS.rhs(),fluidstate_.interfaceSolution());
    interfaceInvOp.finalize();
    timerSolveInterface.stop();

    // set the fluid state for the interface with the correct quantities
    timerSolveInterface.start();
    fluidstate_.finalizeInterfaceQuantities();
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
    std::cout<<"Assemble bulk operators time: "<<timerAssembleBulk.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Assemble interface operator time: "<<timerAssembleInterface.elapsed()<<" seconds."<<std::endl;
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
