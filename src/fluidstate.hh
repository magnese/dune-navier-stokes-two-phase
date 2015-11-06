#ifndef DUNE_FEM_FLUID_STATE_HH
#define DUNE_FEM_FLUID_STATE_HH

// C++ includes
#include <memory>
#include <iostream>
#include <tuple>
#include <algorithm>
#include <string>

// dune includes
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
  OutputParameters(const std::string& prefix):
    prefix_(prefix),startcounter_(0),startcall_(0),startsavetime_(0.0)
  {}

  inline std::string prefix() const
  {
    return prefix_;
  }
  inline int startcounter() const
  {
    return startcounter_;
  }
  inline int startcall() const
  {
    return startcall_;
  }
  inline double startsavetime() const
  {
    return startsavetime_;
  }
  template<typename DOT>
  inline void saveState(const DOT& dataOutput) const
  {
    startcounter_=dataOutput.writeStep();
    startcall_=dataOutput.writeCalls();
    startsavetime_=dataOutput.saveTime();
  }

  private:
  const std::string prefix_;
  mutable int startcounter_;
  mutable int startcall_;
  mutable double startsavetime_;
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
  typedef typename Traits::CurvatureDiscreteSpaceType CurvatureDiscreteSpaceType;
  typedef typename Traits::DisplacementDiscreteSpaceType DisplacementDiscreteSpaceType;
  typedef typename Traits::InterfaceDiscreteSpaceType InterfaceDiscreteSpaceType;
  typedef typename Traits::BulkDisplacementDiscreteSpaceType BulkDisplacementDiscreteSpaceType;
  // define discrete functions
  typedef typename Traits::VelocityDiscreteFunctionType VelocityDiscreteFunctionType;
  typedef typename Traits::PressureDiscreteFunctionType PressureDiscreteFunctionType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef typename Traits::PressureAdditionalDiscreteFunctionType PressureAdditionalDiscreteFunctionType;
  #endif
  typedef typename Traits::PressureDumpDiscreteFunctionType PressureDumpDiscreteFunctionType;
  typedef typename Traits::BulkDiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename Traits::CurvatureDiscreteFunctionType CurvatureDiscreteFunctionType;
  typedef typename Traits::DisplacementDiscreteFunctionType DisplacementDiscreteFunctionType;
  typedef typename Traits::InterfaceDiscreteFunctionType InterfaceDiscreteFunctionType;
  typedef typename Traits::BulkDisplacementDiscreteFunctionType BulkDisplacementDiscreteFunctionType;
  // define data outputs
  typedef std::tuple<const VelocityDiscreteFunctionType*,const PressureDumpDiscreteFunctionType*> BulkTupleType;
  typedef std::tuple<const CurvatureDiscreteFunctionType*,const DisplacementDiscreteFunctionType*> InterfaceTupleType;
  typedef DataOutput<BulkGridType,BulkTupleType> BulkDataOutputType;
  typedef DataOutput<InterfaceGridType,InterfaceTupleType> InterfaceDataOutputType;

  // constructor
  explicit FluidState(int argc,char** argv,const GmshAlgorithmType& algorithm,const bool& verbosity,
                      const bool& checkEntityWithNoVerticesInDomain):
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
      curvaturespace_=other.curvaturespace_;
      displacementspace_=other.displacementspace_;
      interfacespace_=other.interfacespace_;
      bulkdisplacementspace_=other.bulkdisplacementspace_;
      velocity_=other.velocity_;
      pressure_=other.pressure_;
      #if PRESSURE_SPACE_TYPE == 2
      pressureadditional_=other.pressureadditional_;
      pressuredump_=other.pressuredump_;
      #endif
      bulk_=other.bulk_;
      curvature_=other.curvature_;
      displacement_=other.displacement_;
      interface_=other.interface_;
      bulktuple_=other.bulktuple_;
      interfacetuple_=other.interfacetuple_;
      bulkdisplacement_=other.bulkdisplacement_;
      bulkoutput_=other.bulkoutput_;
      interfaceoutput_=other.interfaceoutput_;
      sequence_=other.sequence_;
      bulkoutputparameters_=other.bulkoutputparameters_;
      interfaceoutputparameters_=other.interfaceoutputparameters_;
    }
    return *this;
  }

  // get mesh manager
  inline CoupledMeshManagerType& meshManager()
  {
    return meshmanager_;
  }
  inline const CoupledMeshManagerType& meshManager() const
  {
    return meshmanager_;
  }

  // get grid
  inline BulkGridType& bulkGrid()
  {
    return meshmanager_.bulkGrid();
  }
  inline InterfaceGridType& interfaceGrid()
  {
    return meshmanager_.interfaceGrid();
  }

  // get grid parts
  inline BulkGridPartType& bulkGridPart()
  {
    return *bulkgridpart_;
  }
  inline InterfaceGridPartType& interfaceGridPart()
  {
    return *interfacegridpart_;
  }

  // get spaces
  inline const VelocityDiscreteSpaceType& velocitySpace() const
  {
    return *velocityspace_;
  }
  inline const PressureDiscreteSpaceType& pressureSpace() const
  {
    return *pressurespace_;
  }
  #if PRESSURE_SPACE_TYPE != 2
  inline const PressureDumpDiscreteSpaceType& pressureDumpSpace() const
  {
    return *pressurespace_;
  }
  #else
  inline const PressureAdditionalDiscreteSpaceType& pressureAdditionalSpace() const
  {
    return *pressureadditionalspace_;
  }
  inline const PressureDumpDiscreteSpaceType& pressureDumpSpace() const
  {
    return *pressuredumpspace_;
  }
  #endif
  inline const BulkDiscreteSpaceType& bulkSpace() const
  {
    return *bulkspace_;
  }
  inline const CurvatureDiscreteSpaceType& curvatureSpace() const
  {
    return *curvaturespace_;
  }
  inline const DisplacementDiscreteSpaceType& displacementSpace() const
  {
    return *displacementspace_;
  }
  inline const InterfaceDiscreteSpaceType& interfaceSpace() const
  {
    return *interfacespace_;
  }
  inline const BulkDisplacementDiscreteSpaceType& bulkDisplacementSpace() const
  {
    return *bulkdisplacementspace_;
  }

  // get discrete functions
  inline VelocityDiscreteFunctionType& velocity()
  {
    return *velocity_;
  }
  inline const VelocityDiscreteFunctionType& velocity() const
  {
    return *velocity_;
  }
  inline PressureDiscreteFunctionType& pressure()
  {
    return *pressure_;
  }
  inline const PressureDiscreteFunctionType& pressure() const
  {
    return *pressure_;
  }
  #if PRESSURE_SPACE_TYPE != 2
  inline PressureDumpDiscreteFunctionType& pressureDump()
  {
    return *pressure_;
  }
  inline const PressureDumpDiscreteFunctionType& pressureDump() const
  {
    return *pressure_;
  }
  #else
  inline PressureAdditionalDiscreteFunctionType& pressureAdditional()
  {
    return *pressureadditional_;
  }
  inline const PressureAdditionalDiscreteFunctionType& pressureAdditional() const
  {
    return *pressureadditional_;
  }
  inline PressureDumpDiscreteFunctionType& pressureDump()
  {
    return *pressuredump_;
  }
  inline const PressureDumpDiscreteFunctionType& pressureDump() const
  {
    return *pressuredump_;
  }
  #endif
  inline BulkDiscreteFunctionType& bulkSolution()
  {
    return *bulk_;
  }
  inline const BulkDiscreteFunctionType& bulkSolution() const
  {
    return *bulk_;
  }
  inline CurvatureDiscreteFunctionType& curvature()
  {
    return *curvature_;
  }
  inline const CurvatureDiscreteFunctionType& curvature() const
  {
    return *curvature_;
  }
  inline DisplacementDiscreteFunctionType& displacement()
  {
    return *displacement_;
  }
  inline const DisplacementDiscreteFunctionType& displacement() const
  {
    return *displacement_;
  }
  inline InterfaceDiscreteFunctionType& interfaceSolution()
  {
    return *interface_;
  }
  inline const InterfaceDiscreteFunctionType& interfaceSolution() const
  {
    return *interface_;
  }
  inline BulkDisplacementDiscreteFunctionType& bulkDisplacement()
  {
    return *bulkdisplacement_;
  }
  inline const BulkDisplacementDiscreteFunctionType& bulkDisplacement() const
  {
    return *bulkdisplacement_;
  }

  // rebuild all quantities if the mesh is changed
  bool update()
  {
    bool meshIsChanged(false);
    // check if the mesh is changed
    if(sequence_!=meshmanager_.sequence())
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
      curvaturespace_=std::make_shared<CurvatureDiscreteSpaceType>(interfaceGridPart());
      displacementspace_=std::make_shared<DisplacementDiscreteSpaceType>(interfaceGridPart());
      interfacespace_=std::make_shared<InterfaceDiscreteSpaceType>(interfaceGridPart());
      bulkdisplacementspace_=std::make_shared<BulkDisplacementDiscreteSpaceType>(bulkGridPart());
      // create discrete functions
      velocity_=std::make_shared<VelocityDiscreteFunctionType>("velocity",velocitySpace());
      pressure_=std::make_shared<PressureDiscreteFunctionType>("pressure",pressureSpace());
      #if PRESSURE_SPACE_TYPE == 2
      pressureadditional_=std::make_shared<PressureAdditionalDiscreteFunctionType>("pressure additional",pressureAdditionalSpace());
      pressuredump_=std::make_shared<PressureDumpDiscreteFunctionType>("pressure",pressureDumpSpace());
      #endif
      bulk_=std::make_shared<BulkDiscreteFunctionType>("bulk solution",bulkSpace());
      curvature_=std::make_shared<CurvatureDiscreteFunctionType>("curvature",curvatureSpace());
      displacement_=std::make_shared<DisplacementDiscreteFunctionType>("displacement",displacementSpace());
      interface_=std::make_shared<InterfaceDiscreteFunctionType>("interface solution",interfaceSpace());
      bulkdisplacement_=std::make_shared<BulkDisplacementDiscreteFunctionType>("bulk displacement",bulkDisplacementSpace());
      // create IO
      bulktuple_=std::make_shared<BulkTupleType>(&velocity(),&pressureDump());
      interfacetuple_=std::make_shared<InterfaceTupleType>(&curvature(),&displacement());
      bulkoutput_=std::make_shared<BulkDataOutputType>(meshmanager_.bulkGrid(),*bulktuple_,bulkoutputparameters_);
      interfaceoutput_=std::make_shared<InterfaceDataOutputType>(meshmanager_.interfaceGrid(),*interfacetuple_,interfaceoutputparameters_);
      // print some info
      printInfo();
      // update sequence number
      sequence_=meshmanager_.sequence();
      meshIsChanged=true;
    }
    return meshIsChanged;
  }

  // dump solutions on file
  template<typename TimeProviderType>
  inline void dumpBulkSolutions(const TimeProviderType& timeProvider) const
  {
    #if PRESSURE_SPACE_TYPE == 2
    for(const auto entity:entities(pressure()))
    {
      auto localPressure(pressure().localFunction(entity));
      auto localPressureAdditional(pressureAdditional().localFunction(entity));
      auto localPressureDump(pressureDump().localFunction(entity));
      for(auto i=0;i!=localPressureDump.size();++i)
        localPressureDump[i]=localPressure[i]+localPressureAdditional[0];
    }
    #endif
    bulkoutput_->write(timeProvider);
    bulkoutputparameters_.saveState(*bulkoutput_);
  }
  template<typename TimeProviderType>
  inline void dumpInterfaceSolutions(const TimeProviderType& timeProvider) const
  {
    interfaceoutput_->write(timeProvider);
    interfaceoutputparameters_.saveState(*interfaceoutput_);
  }

  // print spaces info
  const void printInfo(std::ostream& s=std::cout) const
  {
    s<<"Solving for "<<velocity().size()<<" "<<velocity().name()<<" dofs (P"<<velocitySpace().order()<<")."<<std::endl;
    s<<"Solving for "<<pressure().size()<<" "<<pressure().name() <<" dofs (P"<<pressureSpace().order()<<")."<<std::endl;
    #if PRESSURE_SPACE_TYPE == 2
    s<<"Solving for "<<pressure().size()<<" "<<pressureAdditional().name()<<" dofs (P"<<pressureAdditionalSpace().order()<<")."<<std::endl;
    #endif
    s<<"Solving for "<<curvature().size()<<" "<<curvature().name()<<" dofs (P"<<curvatureSpace().order()<<")."<<std::endl;
    s<<"Solving for "<<displacement().size()<<" "<<displacement().name()<<" dofs (P"<<displacementSpace().order()<<")."<<std::endl;
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

  // split interface solution into curvature and displacement
  void finalizeInterfaceQuantities()
  {
    auto solutionIt(interfaceSolution().dbegin());
    for(auto it=curvature().dbegin();it!=curvature().dend();++it,++solutionIt)
      (*it)=(*solutionIt);
    std::copy(solutionIt,interfaceSolution().dend(),displacement().dbegin());
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
    s<<"CurvatureDiscreteSpace = "<<curvaturespace_.use_count()<<std::endl;
    s<<"DisplacementDiscreteSpace = "<<displacementspace_.use_count()<<std::endl;
    s<<"InterfaceDiscreteSpace = "<<interfacespace_.use_count()<<std::endl;
    s<<"BulkDisplacementDiscreteSpace = "<<bulkdisplacementspace_.use_count()<<std::endl;
    s<<"VelocityDiscreteFunction = "<<velocity_.use_count()<<std::endl;
    s<<"PressureDiscreteFunction = "<<pressure_.use_count()<<std::endl;
    #if PRESSURE_SPACE_TYPE == 2
    s<<"PressureAdditionalDiscreteFunction = "<<pressureadditional_.use_count()<<std::endl;
    s<<"PressureDumpDiscreteFunction = "<<pressuredump_.use_count()<<std::endl;
    #endif
    s<<"BulkDiscreteFunction = "<<bulk_.use_count()<<std::endl;
    s<<"CurvatureDiscreteFunction = "<<curvature_.use_count()<<std::endl;
    s<<"DisplacementDiscreteFunction = "<<displacement_.use_count()<<std::endl;
    s<<"InterfaceDiscreteFunction = "<<interface_.use_count()<<std::endl;
    s<<"BulkTuple = "<<bulktuple_.use_count()<<std::endl;
    s<<"InterfaceTuple = "<<interfacetuple_.use_count()<<std::endl;
    s<<"BulkDisplacementDiscreteFunction = "<<bulkdisplacement_.use_count()<<std::endl;
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
  std::shared_ptr<CurvatureDiscreteSpaceType> curvaturespace_;
  std::shared_ptr<DisplacementDiscreteSpaceType> displacementspace_;
  std::shared_ptr<InterfaceDiscreteSpaceType> interfacespace_;
  std::shared_ptr<BulkDisplacementDiscreteSpaceType> bulkdisplacementspace_;
  std::shared_ptr<VelocityDiscreteFunctionType> velocity_;
  std::shared_ptr<PressureDiscreteFunctionType> pressure_;
  #if PRESSURE_SPACE_TYPE == 2
  std::shared_ptr<PressureAdditionalDiscreteFunctionType> pressureadditional_;
  std::shared_ptr<PressureDumpDiscreteFunctionType> pressuredump_;
  #endif
  std::shared_ptr<BulkDiscreteFunctionType> bulk_;
  std::shared_ptr<CurvatureDiscreteFunctionType> curvature_;
  std::shared_ptr<DisplacementDiscreteFunctionType> displacement_;
  std::shared_ptr<InterfaceDiscreteFunctionType> interface_;
  std::shared_ptr<BulkTupleType> bulktuple_;
  std::shared_ptr<InterfaceTupleType> interfacetuple_;
  std::shared_ptr<BulkDisplacementDiscreteFunctionType> bulkdisplacement_;
  std::shared_ptr<BulkDataOutputType> bulkoutput_;
  std::shared_ptr<InterfaceDataOutputType> interfaceoutput_;
  unsigned int sequence_;
  mutable OutputParameters bulkoutputparameters_;
  mutable OutputParameters interfaceoutputparameters_;

  // reset pointers to avoid dangling references
  void resetPointers()
  {
    interfaceoutput_.reset();
    bulkoutput_.reset();
    bulkdisplacement_.reset();
    interfacetuple_.reset();
    bulktuple_.reset();
    interface_.reset();
    displacement_.reset();
    curvature_.reset();
    bulk_.reset();
    #if PRESSURE_SPACE_TYPE == 2
    pressuredump_.reset();
    pressureadditional_.reset();
    #endif
    pressure_.reset();
    velocity_.reset();
    bulkdisplacementspace_.reset();
    interfacespace_.reset();
    displacementspace_.reset();
    curvaturespace_.reset();
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
#endif
