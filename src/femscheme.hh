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
#include <dune/fem/solver/timeprovider.hh>

// local includes
#include "femtraits.hh"
#include "fluidstate.hh"

#include "problems.hh"
#include "nulloperator.hh"
#include "nullrhs.hh"
#include "massmatrix.hh"
#include "operatorwrapper.hh"
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

template<class CoupledMeshManagerImp,class Traits=FemTraits<typename CoupledMeshManagerImp::BulkGridType,
                                                            typename CoupledMeshManagerImp::InterfaceGridType>,
                                     class TimeProviderImp=FixedStepTimeProvider<>>
class FemScheme
{
  public:
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef TimeProviderImp TimeProviderType;
  typedef FemScheme<CoupledMeshManagerType,Traits,TimeProviderType> ThisType;

  // define fluid states
  typedef FluidState<CoupledMeshManagerType,Traits> FluidStateType;

  // define grid parts
  typedef typename Traits::BulkGridType BulkGridType;
  typedef typename Traits::BulkGridPartType BulkGridPartType;
  typedef typename Traits::InterfaceGridType InterfaceGridType;
  typedef typename Traits::InterfaceGridPartType InterfaceGridPartType;

  // define spaces
  typedef typename Traits::VelocityDiscreteSpaceType VelocityDiscreteSpaceType;
  typedef typename Traits::PressureDiscreteSpaceType PressureDiscreteSpaceType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef typename Traits::PressureAdditionalDiscreteSpaceType PressureAdditionalDiscreteSpaceType;
  #endif
  typedef typename Traits::PressureDumpDiscreteSpaceType PressureDumpDiscreteSpaceType;
  typedef typename Traits::CurvatureDiscreteSpaceType CurvatureDiscreteSpaceType;
  typedef typename Traits::DisplacementDiscreteSpaceType DisplacementDiscreteSpaceType;
  typedef typename Traits::BulkDiscreteSpaceType BulkDiscreteSpaceType;
  typedef typename Traits::InterfaceDiscreteSpaceType InterfaceDiscreteSpaceType;

  // define discrete functions
  typedef typename Traits::VelocityDiscreteFunctionType VelocityDiscreteFunctionType;
  typedef typename Traits::PressureDiscreteFunctionType PressureDiscreteFunctionType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef typename Traits::PressureAdditionalDiscreteFunctionType PressureAdditionalDiscreteFunctionType;
  #endif
  typedef typename Traits::PressureDumpDiscreteFunctionType PressureDumpDiscreteFunctionType;
  typedef typename Traits::CurvatureDiscreteFunctionType CurvatureDiscreteFunctionType;
  typedef typename Traits::DisplacementDiscreteFunctionType DisplacementDiscreteFunctionType;
  typedef typename Traits::BulkDiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename Traits::InterfaceDiscreteFunctionType InterfaceDiscreteFunctionType;

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

  // define bulk interface mapper
  typedef typename CoupledMeshManagerType::BulkInterfaceGridMapperType BulkInterfaceGridMapperType;

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
  typedef BulkVelocityOperator<VelocityLinearOperatorType,TimeProviderType,ProblemType> VelocityOperatorType;
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
  typedef BulkVelocityRHS<VelocityDiscreteFunctionType,TimeProviderType,ProblemType> VelocityRHSType;
  #if PROBLEM_NUMBER == 0 || PROBLEM_NUMBER == 1 || PROBLEM_NUMBER == 2 || PROBLEM_NUMBER == 5 || PROBLEM_NUMBER == 6 || PROBLEM_NUMBER == 7
  typedef NullRHS<PressureDiscreteFunctionType,TimeProviderType> PressureRHSType;
  #else
  typedef BulkPressureRHS<PressureDiscreteFunctionType,TimeProviderType> PressureRHSType;
  #endif

  // define operators for interface
  typedef InterfaceOperator<InterfaceLinearOperatorType,TimeProviderType> InterfaceOperatorType;
  typedef VelocityCurvatureOperator<VelocityCurvatureLinearOperatorType,BulkInterfaceGridMapperType> VelocityCurvatureOperatorType;
  typedef InterfaceRHS<InterfaceDiscreteFunctionType> InterfaceRHSType;

  // constructor
  explicit FemScheme(CoupledMeshManagerType& meshManager):
    meshmanager_(meshManager),problem_(meshmanager_)
  {}

  FemScheme(const ThisType& )=delete;

  // obtain mapper
  inline const BulkInterfaceGridMapperType& bulkInterfaceGridMapper() const
  {
    return meshmanager_.mapper();
  }

  // obtain problem
  inline ProblemType& problem()
  {
    return problem_;
  }

  // compute interface initial curvature
  void computeInterface(FluidStateType& fluidState,const TimeProviderType& timeProvider)
  {
    // rebuild all quantities if the mesh is changed
    fluidState.update();

    // assemble operator
    InterfaceOperatorType interfaceOp(fluidState.interfaceSpace(),timeProvider);
    interfaceOp.assemble(bulkInterfaceGridMapper());

    //assemble rhs
    InterfaceRHSType interfaceRHS(fluidState.interfaceSpace());
    interfaceRHS.assemble(interfaceOp);
    interfaceRHS.rhs()*=-1.0;

    // solve
    typedef UMFPACKOp<InterfaceDiscreteFunctionType,InterfaceOperatorType> InterfaceInverseOperatorType;
    InterfaceInverseOperatorType interfaceInvOp(interfaceOp);
    interfaceInvOp(interfaceRHS.rhs(),fluidState.interfaceSolution());

    // set the fluid state for the interface with the correct quantities
    fluidState.finalizeInterfaceQuantities();
  }

  // compute solution
  void operator()(FluidStateType& fluidState,const TimeProviderType& timeProvider)
  {
    // rebuild all quantities if the mesh is changed
    fluidState.update();

    // create timers
    Timer timerAssembleBulk(false);
    Timer timerAssembleInterface(false);
    Timer timerSolveBulk(false);
    Timer timerSolveInterface(false);

    // assemble bulk operators and impose BC
    timerAssembleBulk.start();
    VelocityOperatorType velocityOp(fluidState.velocitySpace(),timeProvider,problem_,fluidState.velocity());
    velocityOp.assemble();
    PressureVelocityOperatorType pressureVelocityOp(fluidState.pressureSpace(),fluidState.velocitySpace());
    pressureVelocityOp.assemble();
    #if PRESSURE_SPACE_TYPE == 2
    PressureAdditionalVelocityOperatorType pressureAdditionalVelocityOp(fluidState.pressureAdditionalSpace(),fluidState.velocitySpace());
    pressureAdditionalVelocityOp.assemble();
    problem_.applyBCToOperator(velocityOp,pressureVelocityOp,pressureAdditionalVelocityOp);
    #else
    problem_.applyBCToOperator(velocityOp,pressureVelocityOp);
    #endif
    VelocityPressureOperatorType velocityPressureOp(fluidState.velocitySpace(),fluidState.pressureSpace());
    velocityPressureOp.assemble();
    #if PRESSURE_SPACE_TYPE == 2
    VelocityPressureAdditionalOperatorType velocityPressureAdditionalOp(fluidState.velocitySpace(),fluidState.pressureAdditionalSpace());
    velocityPressureAdditionalOp.assemble();
    #endif
    PressureOperatorType pressureOp(fluidState.pressureSpace(),fluidState.pressureSpace());
    //pressureOp.assemble(); // not needed since it is never used
    CurvatureVelocityOperatorType curvatureVelocityOp(fluidState.curvatureSpace(),fluidState.velocitySpace(),bulkInterfaceGridMapper());
    curvatureVelocityOp.assemble();
    BulkMassMatrixOperatorType bulkMassMatrixOp(fluidState.pressureSpace());
    bulkMassMatrixOp.assemble();
    #if PRESSURE_SPACE_TYPE == 2
    BulkMassMatrixAdditionalOperatorType bulkMassMatrixAdditionalOp(fluidState.pressureAdditionalSpace());
    bulkMassMatrixAdditionalOp.assemble();
    #endif
    timerAssembleBulk.stop();

    // assemble interface operators
    timerAssembleInterface.start();
    InterfaceOperatorType interfaceOp(fluidState.interfaceSpace(),timeProvider);
    interfaceOp.assemble(bulkInterfaceGridMapper());
    VelocityCurvatureOperatorType velocityCurvatureOp(fluidState.velocitySpace(),fluidState.curvatureSpace(),bulkInterfaceGridMapper());
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
    VelocityRHSType velocityRHS(fluidState.velocitySpace(),timeProvider,problem_,fluidState.velocity());
    velocityRHS.assemble(problem_.velocityRHS());
    PressureRHSType pressureRHS(fluidState.pressureSpace(),timeProvider);
    pressureRHS.assemble(problem_.velocityBC());
    timerAssembleBulk.stop();

    // assemble interface RHS
    timerAssembleInterface.start();
    InterfaceRHSType interfaceRHS(fluidState.interfaceSpace());
    interfaceRHS.assemble(interfaceOp);
    timerAssembleInterface.stop();

    // add bulk coupling
    const auto gamma(problem_.gamma());
    if(gamma!=0.0)
    {
      timerAssembleBulk.start();
      InterfaceDiscreteFunctionType interfaceTempFunction("interface temporary function",fluidState.interfaceSpace());
      interfaceInvOp.apply(interfaceRHS.rhs(),interfaceTempFunction);

      CurvatureDiscreteFunctionType curvatureTempFunction("curvature temporary function",fluidState.curvatureSpace());
      std::copy_n(interfaceTempFunction.dbegin(),curvatureTempFunction.size(),curvatureTempFunction.dbegin());

      VelocityDiscreteFunctionType velocityCouplingRHS("velocity coupling RHS",fluidState.velocitySpace());
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
    BulkDiscreteFunctionType bulkRHS("bulk RHS",fluidState.bulkSpace());
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
    typedef StokesPrecond<VelocityOperatorType,PressureVelocityOperatorType,BulkMassMatrixOperatorType> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(velocityOp,pressureVelocityOp,bulkMassMatrixOp);
    #else
    typedef ExtendedOperatorWrapper<CoupledOperatorWrapperType,PressureVelocityOperatorType,VelocityPressureOperatorType,
                                    PressureAdditionalVelocityOperatorType,VelocityPressureAdditionalOperatorType> BulkOperatorWrapperType;
    BulkOperatorWrapperType bulkOp(coupledWrapperOp,pressureVelocityOp,velocityPressureOp,pressureAdditionalVelocityOp,
                                   velocityPressureAdditionalOp);
    typedef ExtendedStokesPrecond<VelocityOperatorType,PressureVelocityOperatorType,BulkMassMatrixOperatorType,
                                  PressureAdditionalVelocityOperatorType,BulkMassMatrixAdditionalOperatorType> BulkPreconditionerType;
    BulkPreconditionerType bulkPreconditioner(velocityOp,pressureVelocityOp,bulkMassMatrixOp,pressureAdditionalVelocityOp,
                                              bulkMassMatrixAdditionalOp);
    #endif

    typedef Dune::RestartedGMResSolver<typename BulkDiscreteFunctionType::DofStorageType> BulkLinearInverseOperatorType;
    InverseOperatorResult returnInfo;
    BulkLinearInverseOperatorType bulkInvOp(bulkOp,bulkPreconditioner,redEps,restart,maxIter,verbosity);
    bulkInvOp.apply(fluidState.bulkSolution().blockVector(),bulkRHS.blockVector(),returnInfo);
    timerSolveBulk.stop();

    // set the fluid state for the bulk with the correct quantities
    fluidState.finalizeBulkQuantities();
    timerSolveBulk.stop();

    // add interface coupling
    timerAssembleInterface.start();
    CurvatureDiscreteFunctionType curvatureCouplingRHS("curvature coupling RHS",fluidState.curvatureSpace());
    velocityCurvatureOp(fluidState.velocity(),curvatureCouplingRHS);
    std::copy(curvatureCouplingRHS.dbegin(),curvatureCouplingRHS.dend(),interfaceRHS.rhs().dbegin());
    interfaceRHS.rhs()*=-1.0;
    timerAssembleInterface.stop();

    // solve interface
    timerSolveInterface.start();
    interfaceInvOp.apply(interfaceRHS.rhs(),fluidState.interfaceSolution());
    interfaceInvOp.finalize();
    timerSolveInterface.stop();

    // set the fluid state for the interface with the correct quantities
    timerSolveInterface.start();
    fluidState.finalizeInterfaceQuantities();
    timerSolveInterface.stop();

    // project pressure solution to the space of mean zero function
    timerSolveBulk.start();
    projectZeroMean(fluidState.pressure(),bulkMassMatrixOp);
    #if PRESSURE_SPACE_TYPE == 2
    projectZeroMean(fluidState.pressureAdditional(),bulkMassMatrixAdditionalOp);
    #endif
    timerSolveBulk.stop();

    // remove displacement if the problem has only 1 phase
    if(gamma==0.0)
      fluidState.displacement().clear();

    // print timers
    std::cout<<"Assemble bulk operators time: "<<timerAssembleBulk.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Assemble interface operator time: "<<timerAssembleInterface.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Solve bulk time: "<<timerSolveBulk.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Solve interface time: "<<timerSolveInterface.elapsed()<<" seconds."<<std::endl;
  }

  void remesh()
  {
    // create timer
    Timer timer(false);
    // remesh
    timer.start();
    auto bulkHostGrid(meshmanager_.remesh());
    // finalize remeshing
    meshmanager_.finalize(bulkHostGrid);
    timer.stop();
    if(bulkHostGrid!=nullptr)
      std::cout<<"Remeshing time: "<<timer.elapsed()<<" seconds."<<std::endl;
  }

  private:
  CoupledMeshManagerType& meshmanager_;
  ProblemType problem_;

  // project df into the space of mean zero function
  template<class DF,class Op>
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
