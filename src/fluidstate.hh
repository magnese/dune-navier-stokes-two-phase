#ifndef DUNE_FEM_FLUIDSTATE_HH
#define DUNE_FEM_FLUIDSTATE_HH

// C++ includes
#include <memory>
#include <iostream>
#include <tuple>
#include <algorithm>
#include <string>

// dune includes
#include <dune/fem/common/tupleforeach.hh>
#include <dune/fem/io/file/dataoutput.hh>

// local includes
#include "femtraits.hh"

namespace Dune
{
namespace Fem
{

class OutputParameters:public DataOutputParameters
{
  public:
  OutputParameters(std::string&& prefix):
    prefix_(prefix),startcounter_(0),startcall_(0),startsavetime_(0.0)
  {}

  std::string prefix() const
  {
    return prefix_;
  }
  int startcounter() const
  {
    return startcounter_;
  }
  int startcall() const
  {
    return startcall_;
  }
  double startsavetime() const
  {
    return startsavetime_;
  }
  template<typename DOT>
  void saveState(const DOT& dataOutput)
  {
    startcounter_=dataOutput.writeStep();
    startcall_=dataOutput.writeCalls();
    startsavetime_=dataOutput.saveTime();
  }

  private:
  const std::string prefix_;
  int startcounter_;
  int startcall_;
  double startsavetime_;
};

template<typename CoupledMeshManagerImp,typename Traits=FemTraits<typename CoupledMeshManagerImp::BulkGridType,
                                                                  typename CoupledMeshManagerImp::InterfaceGridType>>
class FluidState
{
  public:
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef FluidState<CoupledMeshManagerType,Traits> ThisType;

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
  typedef typename Traits::BulkDiscreteSpaceType BulkDiscreteSpaceType;
  typedef typename Traits::BulkDisplacementDiscreteSpaceType BulkDisplacementDiscreteSpaceType;
  typedef typename Traits::CurvatureDiscreteSpaceType CurvatureDiscreteSpaceType;
  typedef typename Traits::DisplacementDiscreteSpaceType DisplacementDiscreteSpaceType;
  typedef typename Traits::InterfaceDiscreteSpaceType InterfaceDiscreteSpaceType;
  // define discrete functions
  typedef typename Traits::VelocityDiscreteFunctionType VelocityDiscreteFunctionType;
  typedef typename Traits::PressureDiscreteFunctionType PressureDiscreteFunctionType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef typename Traits::PressureAdditionalDiscreteFunctionType PressureAdditionalDiscreteFunctionType;
  #endif
  typedef typename Traits::PressureDumpDiscreteFunctionType PressureDumpDiscreteFunctionType;
  typedef typename Traits::BulkDiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename Traits::BulkDisplacementDiscreteFunctionType BulkDisplacementDiscreteFunctionType;
  typedef typename Traits::CurvatureDiscreteFunctionType CurvatureDiscreteFunctionType;
  typedef typename Traits::DisplacementDiscreteFunctionType DisplacementDiscreteFunctionType;
  typedef typename Traits::InterfaceDiscreteFunctionType InterfaceDiscreteFunctionType;
  // define data outputs
  typedef std::tuple<VelocityDiscreteFunctionType*,PressureDumpDiscreteFunctionType*> BulkTupleType;
  typedef std::tuple<CurvatureDiscreteFunctionType*> InterfaceTupleType;
  typedef DataOutput<BulkGridType,BulkTupleType> BulkDataOutputType;
  typedef DataOutput<InterfaceGridType,InterfaceTupleType> InterfaceDataOutputType;

  // constructor
  explicit FluidState(int argc,char** argv,const GmshAlgorithmType& algorithm,bool verbosity,bool checkEntityWithNoVerticesInDomain):
    meshmanager_(argc,argv,algorithm,verbosity,checkEntityWithNoVerticesInDomain),sequence_(0),bulkoutputparameters_("bulk-"),
    interfaceoutputparameters_("interface-")
  {}

  // destructor
  ~FluidState()
  {
    resetPointers();
  }

  // copy constructor
  FluidState(const ThisType& )=default;

  // copy assignament
  ThisType& operator=(const ThisType& other)
  {
    if(this!=&other)
    {
      resetPointers();
      meshmanager_=other.meshmanager_;
      bulkgridpart_=other.bulkgridpart_;
      interfacegridpart_=other.interfacegridpart_;
      velocityspace_=other.velocityspace_;
      pressurespace_=other.pressurespace_;
      #if PRESSURE_SPACE_TYPE == 2
      pressureadditionalspace_=other.pressureadditionalspace_;
      pressuredumpspace_=other.pressuredumpspace_;
      #endif
      bulkspace_=other.bulkspace_;
      bulkdisplacementspace_=other.bulkdisplacementspace_;
      interfacespace_=other.interfacespace_;
      velocity_=other.velocity_;
      pressure_=other.pressure_;
      #if PRESSURE_SPACE_TYPE == 2
      pressureadditional_=other.pressureadditional_;
      pressuredump_=other.pressuredump_;
      #endif
      bulk_=other.bulk_;
      bulkdisplacement_=other.bulkdisplacement_;
      interface_=other.interface_;
      bulktuple_=other.bulktuple_;
      interfacetuple_=other.interfacetuple_;
      bulkoutput_=other.bulkoutput_;
      interfaceoutput_=other.interfaceoutput_;
      sequence_=other.sequence_;
      bulkoutputparameters_=other.bulkoutputparameters_;
      interfaceoutputparameters_=other.interfaceoutputparameters_;
    }
    return *this;
  }

  // get mesh manager
  CoupledMeshManagerType& meshManager()
  {
    return meshmanager_;
  }
  const CoupledMeshManagerType& meshManager() const
  {
    return meshmanager_;
  }

  // get grid
  BulkGridType& bulkGrid()
  {
    return meshmanager_.bulkGrid();
  }
  InterfaceGridType& interfaceGrid()
  {
    return meshmanager_.interfaceGrid();
  }

  // get grid parts
  BulkGridPartType& bulkGridPart()
  {
    return *bulkgridpart_;
  }
  InterfaceGridPartType& interfaceGridPart()
  {
    return *interfacegridpart_;
  }

  // get spaces
  const VelocityDiscreteSpaceType& velocitySpace() const
  {
    return *velocityspace_;
  }
  const PressureDiscreteSpaceType& pressureSpace() const
  {
    return *pressurespace_;
  }
  #if PRESSURE_SPACE_TYPE != 2
  const PressureDumpDiscreteSpaceType& pressureDumpSpace() const
  {
    return *pressurespace_;
  }
  #else
  const PressureAdditionalDiscreteSpaceType& pressureAdditionalSpace() const
  {
    return *pressureadditionalspace_;
  }
  const PressureDumpDiscreteSpaceType& pressureDumpSpace() const
  {
    return *pressuredumpspace_;
  }
  #endif
  const BulkDiscreteSpaceType& bulkSpace() const
  {
    return *bulkspace_;
  }
  const BulkDisplacementDiscreteSpaceType& bulkDisplacementSpace() const
  {
    return *bulkdisplacementspace_;
  }
  const CurvatureDiscreteSpaceType& curvatureSpace() const
  {
    return interfacespace_->template subDiscreteFunctionSpace<0>();
  }
  const DisplacementDiscreteSpaceType& displacementSpace() const
  {
    return interfacespace_->template subDiscreteFunctionSpace<1>();
  }
  const InterfaceDiscreteSpaceType& interfaceSpace() const
  {
    return *interfacespace_;
  }

  // get discrete functions
  VelocityDiscreteFunctionType& velocity()
  {
    return *velocity_;
  }
  const VelocityDiscreteFunctionType& velocity() const
  {
    return *velocity_;
  }
  PressureDiscreteFunctionType& pressure()
  {
    return *pressure_;
  }
  const PressureDiscreteFunctionType& pressure() const
  {
    return *pressure_;
  }
  #if PRESSURE_SPACE_TYPE != 2
  PressureDumpDiscreteFunctionType& pressureDump()
  {
    return *pressure_;
  }
  const PressureDumpDiscreteFunctionType& pressureDump() const
  {
    return *pressure_;
  }
  #else
  PressureAdditionalDiscreteFunctionType& pressureAdditional()
  {
    return *pressureadditional_;
  }
  const PressureAdditionalDiscreteFunctionType& pressureAdditional() const
  {
    return *pressureadditional_;
  }
  PressureDumpDiscreteFunctionType& pressureDump()
  {
    return *pressuredump_;
  }
  const PressureDumpDiscreteFunctionType& pressureDump() const
  {
    return *pressuredump_;
  }
  #endif
  BulkDiscreteFunctionType& bulkSolution()
  {
    return *bulk_;
  }
  const BulkDiscreteFunctionType& bulkSolution() const
  {
    return *bulk_;
  }
  BulkDisplacementDiscreteFunctionType& bulkDisplacement()
  {
    return *bulkdisplacement_;
  }
  const BulkDisplacementDiscreteFunctionType& bulkDisplacement() const
  {
    return *bulkdisplacement_;
  }
  CurvatureDiscreteFunctionType& curvature()
  {
    return interface_->template subDiscreteFunction<0>();
  }
  const CurvatureDiscreteFunctionType& curvature() const
  {
    return interface_->template subDiscreteFunction<0>();
  }
  DisplacementDiscreteFunctionType& displacement()
  {
    return interface_->template subDiscreteFunction<1>();
  }
  const DisplacementDiscreteFunctionType& displacement() const
  {
    return interface_->template subDiscreteFunction<1>();
  }
  InterfaceDiscreteFunctionType& interfaceSolution()
  {
    return *interface_;
  }
  const InterfaceDiscreteFunctionType& interfaceSolution() const
  {
    return *interface_;
  }

  // rebuild all quantities
  void init()
  {
    // create grid parts
    resetPointers();
    bulkgridpart_=std::make_shared<BulkGridPartType>(bulkGrid());
    interfacegridpart_=std::make_shared<InterfaceGridPartType>(interfaceGrid());
    // create spaces
    velocityspace_=std::make_shared<VelocityDiscreteSpaceType>(bulkGridPart());
    pressurespace_=std::make_shared<PressureDiscreteSpaceType>(bulkGridPart());
    #if PRESSURE_SPACE_TYPE == 2
    pressureadditionalspace_=std::make_shared<PressureAdditionalDiscreteSpaceType>(bulkGridPart());
    pressuredumpspace_=std::make_shared<PressureDumpDiscreteSpaceType>(bulkGridPart());
    #endif
    bulkspace_=std::make_shared<BulkDiscreteSpaceType>(bulkGridPart());
    bulkdisplacementspace_=std::make_shared<BulkDisplacementDiscreteSpaceType>(bulkGridPart());
    interfacespace_=std::make_shared<InterfaceDiscreteSpaceType>(interfaceGridPart());
    // create discrete functions
    velocity_=std::make_shared<VelocityDiscreteFunctionType>("velocity",velocitySpace());
    pressure_=std::make_shared<PressureDiscreteFunctionType>("pressure",pressureSpace());
    #if PRESSURE_SPACE_TYPE == 2
    pressureadditional_=std::make_shared<PressureAdditionalDiscreteFunctionType>("pressure additional",pressureAdditionalSpace());
    pressuredump_=std::make_shared<PressureDumpDiscreteFunctionType>("pressure",pressureDumpSpace());
    #endif
    bulk_=std::make_shared<BulkDiscreteFunctionType>("bulk solution",bulkSpace());
    bulkdisplacement_=std::make_shared<BulkDisplacementDiscreteFunctionType>("bulk displacement",bulkDisplacementSpace());
    interface_=std::make_shared<InterfaceDiscreteFunctionType>("interface solution",interfaceSpace());
    // create IO
    bulktuple_=std::make_tuple(&velocity(),&pressureDump());
    interfacetuple_=std::make_tuple(&curvature());
    bulkoutput_=std::make_shared<BulkDataOutputType>(meshmanager_.bulkGrid(),bulktuple_,bulkoutputparameters_);
    interfaceoutput_=std::make_shared<InterfaceDataOutputType>(meshmanager_.interfaceGrid(),interfacetuple_,interfaceoutputparameters_);
    // update sequence number
    sequence_=meshmanager_.sequence();
  }

  // rebuild all quantities if the mesh is changed
  bool update()
  {
    bool meshIsChanged(false);
    // check if the mesh is changed
    if(sequence_!=meshmanager_.sequence())
    {
      init();
      meshIsChanged=true;
    }
    return meshIsChanged;
  }

  // dump solutions on file
  template<typename TimeProviderType>
  void dumpBulkSolutions(const TimeProviderType& timeProvider)
  {
    #if PRESSURE_SPACE_TYPE == 2
    for(const auto& entity:entities(pressure()))
    {
      auto localPressure(pressure().localFunction(entity));
      auto localPressureAdditional(pressureAdditional().localFunction(entity));
      auto localPressureDump(pressureDump().localFunction(entity));
      for(auto i=decltype(localPressureDump.size()){0};i!=localPressureDump.size();++i)
        localPressureDump[i]=localPressure[i]+localPressureAdditional[0];
    }
    #endif
    bulkoutput_->write(timeProvider);
    bulkoutputparameters_.saveState(*bulkoutput_);
  }
  template<typename TimeProviderType>
  void dumpInterfaceSolutions(const TimeProviderType& timeProvider)
  {
    interfaceoutput_->write(timeProvider);
    interfaceoutputparameters_.saveState(*interfaceoutput_);
  }

  // print bulk spaces info
  const void printBulkInfo(std::ostream& s=std::cout) const
  {
    s<<" P"<<velocitySpace().order()<<" "<<velocity().name()<<" -> "<<velocity().size()<<" DOFs ";
    s<<" P"<<pressureSpace().order()<<" "<<pressure().name()<<" -> "<<pressure().size()<<" DOFs ";
    #if PRESSURE_SPACE_TYPE == 2
    s<<" P"<<pressureAdditionalSpace().order()<<" "<<pressureAdditional().name()<<" -> "<<pressureAdditional().size()<<" DOFs ";
    #endif
  }

  // print interface spaces info
  const void printInterfaceInfo(std::ostream& s=std::cout) const
  {
    for_each(*interface_,
      [&s](const auto& block,auto I){s<<" P"<<block.space().order()<<" "<<block.name()<<" -> "<<block.size()<<" DOFs ";});
  }

  // split bulk solution into velocity and pressure
  void finalizeBulkQuantities()
  {
    auto solutionIt(bulkSolution().dbegin());
    for(auto it=velocity().dbegin();it!=velocity().dend();++it,++solutionIt)
      (*it)=(*solutionIt);
    #if PRESSURE_SPACE_TYPE != 2
    std::copy(solutionIt,bulkSolution().dend(),pressure().dbegin());
    #else
    for(auto it=pressure().dbegin();it!=pressure().dend();++it,++solutionIt)
      (*it)=(*solutionIt);
    std::copy(solutionIt,bulkSolution().dend(),pressureAdditional().dbegin());
    #endif
  }

  // check pointers status
  void checkPointersStatus(const std::string& str="FluidState",std::ostream& s=std::cout) const
  {
    s<<std::endl;
    s<<str<<" number of pointers for each object (sequence = "<<sequence_<<") :"<<std::endl;
    s<<"BulkGridPart = "<<bulkgridpart_.use_count()<<std::endl;
    s<<"InterfaceGridPart = "<<interfacegridpart_.use_count()<<std::endl;
    s<<"VelocityDiscreteSpace = "<<velocityspace_.use_count()<<std::endl;
    s<<"PressureDiscreteSpace = "<<pressurespace_.use_count()<<std::endl;
    #if PRESSURE_SPACE_TYPE == 2
    s<<"PressureAdditionalDiscreteSpace = "<<pressureadditionalspace_.use_count()<<std::endl;
    s<<"PressureDumpDiscreteSpace = "<<pressuredumpspace_.use_count()<<std::endl;
    #endif
    s<<"BulkDiscreteSpace = "<<bulkspace_.use_count()<<std::endl;
    s<<"BulkDisplacementDiscreteSpace = "<<bulkdisplacementspace_.use_count()<<std::endl;
    s<<"InterfaceDiscreteSpace = "<<interfacespace_.use_count()<<std::endl;
    s<<"VelocityDiscreteFunction = "<<velocity_.use_count()<<std::endl;
    s<<"PressureDiscreteFunction = "<<pressure_.use_count()<<std::endl;
    #if PRESSURE_SPACE_TYPE == 2
    s<<"PressureAdditionalDiscreteFunction = "<<pressureadditional_.use_count()<<std::endl;
    s<<"PressureDumpDiscreteFunction = "<<pressuredump_.use_count()<<std::endl;
    #endif
    s<<"BulkDiscreteFunction = "<<bulk_.use_count()<<std::endl;
    s<<"BulkDisplacementDiscreteFunction = "<<bulkdisplacement_.use_count()<<std::endl;
    s<<"InterfaceDiscreteFunction = "<<interface_.use_count()<<std::endl;
    s<<"BulkDataOutput = "<<bulkoutput_.use_count()<<std::endl;
    s<<"InterfaceDataOutput = "<<interfaceoutput_.use_count()<<std::endl;
    s<<std::endl;
  }

  private:
  CoupledMeshManagerType meshmanager_;
  std::shared_ptr<BulkGridPartType> bulkgridpart_;
  std::shared_ptr<InterfaceGridPartType> interfacegridpart_;
  std::shared_ptr<VelocityDiscreteSpaceType> velocityspace_;
  std::shared_ptr<PressureDiscreteSpaceType> pressurespace_;
  #if PRESSURE_SPACE_TYPE == 2
  std::shared_ptr<PressureAdditionalDiscreteSpaceType> pressureadditionalspace_;
  std::shared_ptr<PressureDumpDiscreteSpaceType> pressuredumpspace_;
  #endif
  std::shared_ptr<BulkDiscreteSpaceType> bulkspace_;
  std::shared_ptr<BulkDisplacementDiscreteSpaceType> bulkdisplacementspace_;
  std::shared_ptr<InterfaceDiscreteSpaceType> interfacespace_;
  std::shared_ptr<VelocityDiscreteFunctionType> velocity_;
  std::shared_ptr<PressureDiscreteFunctionType> pressure_;
  #if PRESSURE_SPACE_TYPE == 2
  std::shared_ptr<PressureAdditionalDiscreteFunctionType> pressureadditional_;
  std::shared_ptr<PressureDumpDiscreteFunctionType> pressuredump_;
  #endif
  std::shared_ptr<BulkDiscreteFunctionType> bulk_;
  std::shared_ptr<BulkDisplacementDiscreteFunctionType> bulkdisplacement_;
  std::shared_ptr<InterfaceDiscreteFunctionType> interface_;
  InterfaceTupleType interfacetuple_;
  BulkTupleType bulktuple_;
  std::shared_ptr<BulkDataOutputType> bulkoutput_;
  std::shared_ptr<InterfaceDataOutputType> interfaceoutput_;
  unsigned int sequence_;
  OutputParameters bulkoutputparameters_;
  OutputParameters interfaceoutputparameters_;

  // reset pointers to avoid dangling references
  void resetPointers()
  {
    interfaceoutput_.reset();
    bulkoutput_.reset();
    interface_.reset();
    bulkdisplacement_.reset();
    bulk_.reset();
    #if PRESSURE_SPACE_TYPE == 2
    pressuredump_.reset();
    pressureadditional_.reset();
    #endif
    pressure_.reset();
    velocity_.reset();
    interfacespace_.reset();
    bulkdisplacementspace_.reset();
    bulkspace_.reset();
    #if PRESSURE_SPACE_TYPE == 2
    pressuredumpspace_.reset();
    pressureadditionalspace_.reset();
    #endif
    pressurespace_.reset();
    velocityspace_.reset();
    interfacegridpart_.reset();
    bulkgridpart_.reset();
  }
};

}
}

#endif // DUNE_FEM_FLUIDSTATE_HH
