#define GRIDDIM ALBERTA_DIM
#define WORLDDIM ALBERTA_DIM

//#define ALUGRID_SIMPLEX
#define ALBERTAGRID

#define PRESSURE_SPACE_TYPE 0 // 0 P0, 1 P1, 2 P1+P0

#define PRECONDITIONER_TYPE 1 // 0 blocks UMFPACK, 1 full UMFPACK, 2 full SPQR, 3 full LDL

#define REMESH_TYPE 1 // 0 none, 1 uniform, 2 fixed, 3 adaptive

#define INTERPOLATION_TYPE 2 // 0 linear search, 1 barycentric coordinate search, 2 sorted entities + barycentric search

#define PROBLEM_NUMBER 2 // 0 StokesTest1, 1 StokesTest2, 2 StationaryBubble, 3 ExpandingBubble, 4 ShearFlow
                         // 5 StationaryNavierStokes, 6 NavierStokes2D, 7 RisingBubble, 8 NavierStokesTest1, 9 NavierStokesTest2

#if PRECONDITIONER_TYPE == 3
#define USE_SYMMETRIC_DIRICHLET 1
#else
#define USE_SYMMETRIC_DIRICHLET 0
#endif

#include <string>
#include <vector>
#include <memory>
#include <algorithm>

#include "config.h"
#include <dune/common/timer.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/parameter.hh>

#include "coupledmeshmanager.hh"
#include "fluidstate.hh"
#include "femscheme.hh"
#include "meshsmoothing.hh"
#include "compute.hh"

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
    Dune::Fem::Parameter::append((argc<2)?("/home/ma/m/ma2413/dune-repo/dune-navier-stokes-two-phase/src/parameter"):(argv[1]));
    #if REMESH_TYPE != 0
    GmshInitialize(argc,argv);
    #endif

    // create coupled grids and fluid state
    typedef Dune::GridSelector::GridType BulkHostGridType;
    constexpr unsigned int worlddim(BulkHostGridType::dimensionworld);
    constexpr unsigned int bulkGriddim(BulkHostGridType::dimension);
    typedef Dune::AlbertaGrid<bulkGriddim-1,worlddim> InterfaceHostGridType;
    #if REMESH_TYPE == 0
    typedef Dune::Fem::CoupledMeshManager<BulkHostGridType,InterfaceHostGridType,false> CoupledMeshManagerType;
    #elif REMESH_TYPE == 1
    typedef Dune::Fem::CoupledMeshManager<BulkHostGridType,InterfaceHostGridType,true, Dune::UniformCharlength> CoupledMeshManagerType;
    #elif REMESH_TYPE == 2
    typedef Dune::Fem::CoupledMeshManager<BulkHostGridType,InterfaceHostGridType,true, Dune::FixedCharlength> CoupledMeshManagerType;
    #else
    typedef Dune::Fem::CoupledMeshManager<BulkHostGridType,InterfaceHostGridType,true, Dune::AdaptiveCharlength> CoupledMeshManagerType;
    #endif
    #if PRESSURE_SPACE_TYPE == 0
    constexpr bool checkEntityWithNoVerticesInDomain(false);
    #else
    constexpr bool checkEntityWithNoVerticesInDomain(true);
    #endif
    typedef Dune::Fem::FluidState<CoupledMeshManagerType> FluidStateType;
    FluidStateType fluidState(argc,argv,Dune::automatic,false,checkEntityWithNoVerticesInDomain);
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
    std::cout<<std::endl<<"Total running time: "<<timer.elapsed()<<" seconds."<<std::endl;

    return 0;
  }

  catch(std::exception& e)
  {
    throw;
  }

  catch(...)
  {
    std::cerr<<"Unknown exception thrown!"<<std::endl;
    exit(1);
  }
}
