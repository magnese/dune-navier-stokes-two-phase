#define PRESSURE_SPACE_TYPE 0 // 0 P0, 1 P1, 2 P1+P0, 3 P1-DG

#define PRECONDITIONER_TYPE 1 // 0 blocks UMFPACK, 1 full UMFPACK, 2 full SPQR

#define REMESH_TYPE 1 // 0 none, 1 uniform, 2 fixed, 3 adaptive

#define REMESH_CRITERION 0 // 0 volume, 1 angle

#define INTERPOLATION_TYPE 2 // 0 linear search, 1 barycentric coordinate search, 2 sorted entities + barycentric search

#define PROBLEM_NUMBER 7 // 0 StokesTest1, 1 StokesTest2, 2 StationaryBubble, 3 ExpandingBubble, 4 ShearFlow
                         // 5 StationaryNavierStokes, 6 NavierStokes2D, 7 RisingBubble, 8 NavierStokesExpandingBubble1,
                         // 9 NavierStokesExpandingBubble2, 10 NavierStokesExpandingBubble3

#define USE_SYMMETRIC_LAPLACIAN_TERM 1
#define USE_ANTISYMMETRIC_CONVECTIVE_TERM 0
#define USE_SYMMETRIC_DIRICHLET 0

#if PRESSURE_SPACE_TYPE == 0 || PRESSURE_SPACE_TYPE == 1
#define USE_EXTENDED_PRESSURE_SPACE 0
#elif PRESSURE_SPACE_TYPE == 2 || PRESSURE_SPACE_TYPE == 3
#define USE_EXTENDED_PRESSURE_SPACE 1
#endif

#define STRINGIZE_(x) #x
#define STRINGIZE(x) STRINGIZE_(x)

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "config.h"
#include <dune/common/timer.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/grid/albertagrid.hh>

#include "compute.hh"
#include "coupledmeshmanager.hh"
#include "femscheme.hh"
#include "fluidstate.hh"
#include "meshsmoothing.hh"

int main(int argc,char** argv)
{
  try
  {
    // start timer
    Dune::Timer timer(false);
    timer.start();

    // init
    Dune::Fem::MPIManager::initialize(argc,argv);
    Dune::Fem::Parameter::append(argc,argv);
    Dune::Fem::Parameter::append(argc<2?(static_cast<std::string>(SOURCEDIR)+"/parameter"):argv[1]);

    // add preprocessor variables to the parameters in order to have them dumped in the output
    Dune::Fem::Parameter::append("PreconditionerType",STRINGIZE(PRECONDITIONER_TYPE));
    Dune::Fem::Parameter::append("RemeshCriterion",STRINGIZE(REMESH_CRITERION));
    Dune::Fem::Parameter::append("InterpolationType",STRINGIZE(INTERPOLATION_TYPE));
    Dune::Fem::Parameter::append("UseSymmetricLaplacianTerm",STRINGIZE(USE_SYMMETRIC_LAPLACIAN_TERM));
    Dune::Fem::Parameter::append("UseAntisymmetricConvectiveTerm",STRINGIZE(USE_ANTISYMMETRIC_CONVECTIVE_TERM));
    Dune::Fem::Parameter::append("UseSymmetricDirichlet",STRINGIZE(USE_SYMMETRIC_DIRICHLET));

    #if REMESH_TYPE == 1 || REMESH_TYPE == 2 || REMESH_TYPE == 3
    GmshInitialize(argc,argv);
    #endif

    // create coupled grids and fluid state
    constexpr unsigned int worlddim(WORLDDIM);
    constexpr unsigned int bulkGriddim(GRIDDIM);
    typedef Dune::AlbertaGrid<bulkGriddim,worlddim> BulkHostGridType;
    typedef Dune::AlbertaGrid<bulkGriddim-1,worlddim> InterfaceHostGridType;
    #if REMESH_TYPE == 0
    constexpr bool UseCompoundManager(false);
    typedef Dune::UniformCharlength CharlengthPolicyType;
    #elif REMESH_TYPE == 1
    constexpr bool UseCompoundManager(true);
    typedef Dune::UniformCharlength CharlengthPolicyType;
    #elif REMESH_TYPE == 2
    constexpr bool UseCompoundManager(true);
    typedef Dune::FixedCharlength CharlengthPolicyType;
    #elif REMESH_TYPE == 3
    constexpr bool UseCompoundManager(true);
    typedef Dune::AdaptiveCharlength CharlengthPolicyType;
    #endif
    #if REMESH_CRITERION == 0
    typedef Dune::Fem::RemeshingVolumeCriterion RemeshingCriterionType;
    #elif REMESH_CRITERION == 1
    typedef Dune::Fem::RemeshingAngleCriterion RemeshingCriterionType;
    #endif
    typedef Dune::Fem::CoupledMeshManager<BulkHostGridType,InterfaceHostGridType,UseCompoundManager,CharlengthPolicyType,
      RemeshingCriterionType> CoupledMeshManagerType;
    #if PRESSURE_SPACE_TYPE == 0
    constexpr bool checkEntityWithNoVerticesInDomain(false);
    #else
    constexpr bool checkEntityWithNoVerticesInDomain(true);
    #endif
    typedef Dune::Fem::FluidState<CoupledMeshManagerType> FluidStateType;
    FluidStateType fluidState(argc,argv,checkEntityWithNoVerticesInDomain);
    fluidState.meshManager().printInfo();

    // create mesh smoothing
    typedef Dune::Fem::MeshSmoothing<FluidStateType> MeshSmoothingType;
    MeshSmoothingType smoothing(fluidState);
    smoothing.printInfo();

    // compute solution
    typedef Dune::Fem::FemScheme<FluidStateType> FemSchemeType;
    FemSchemeType femScheme(fluidState);
    femScheme.problem().printInfo();
    std::vector<double> errors;
    Dune::Fem::compute<FemSchemeType,MeshSmoothingType>(femScheme,smoothing,errors);

    // output total running time
    timer.stop();
    std::cout<<"\nTotal running time: "<<timer.elapsed()<<" seconds.\n";

    // output parameters
    std::cout<<"\nParameters used:\n";
    Dune::Fem::Parameter::write(std::cout);

    return 0;
  }

  catch(std::exception& e)
  {
    throw;
  }

  catch(...)
  {
    std::cerr<<"Unknown exception thrown!\n";
    exit(1);
  }
}
