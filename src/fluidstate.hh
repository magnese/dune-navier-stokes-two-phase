#ifndef DUNE_FEM_FLUIDSTATE_HH
#define DUNE_FEM_FLUIDSTATE_HH

#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/fem/io/file/dataoutput.hh>

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

template<typename CoupledMeshManagerImp,typename Traits=FemTraits<CoupledMeshManagerImp>>
class FluidState
{
  public:
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef FluidState<CoupledMeshManagerType,Traits> ThisType;

  // define grid parts
  typedef typename Traits::BulkGridType BulkGridType;
  typedef typename Traits::BulkGridPartType BulkGridPartType;
  typedef typename Traits::BulkInnerGridPartType BulkInnerGridPartType;
  typedef typename Traits::BulkOuterGridPartType BulkOuterGridPartType;
  typedef typename Traits::InterfaceGridType InterfaceGridType;
  typedef typename Traits::InterfaceGridPartType InterfaceGridPartType;
  // define spaces
  typedef typename Traits::VelocityDiscreteSpaceType VelocityDiscreteSpaceType;
  typedef typename Traits::PressureDiscreteSpaceType PressureDiscreteSpaceType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef typename Traits::PressureAdditionalDiscreteSpaceType PressureAdditionalDiscreteSpaceType;
  #endif
  typedef typename Traits::PressureDumpDiscreteSpaceType PressureDumpDiscreteSpaceType;
  typedef typename Traits::BulkDisplacementDiscreteSpaceType BulkDisplacementDiscreteSpaceType;
  typedef typename Traits::PhysicalCoefficientDiscreteSpaceType PhysicalCoefficientDiscreteSpaceType;
  typedef typename Traits::CurvatureDiscreteSpaceType CurvatureDiscreteSpaceType;
  typedef typename Traits::DisplacementDiscreteSpaceType DisplacementDiscreteSpaceType;
  // define discrete functions
  typedef typename Traits::VelocityDiscreteFunctionType VelocityDiscreteFunctionType;
  typedef typename Traits::PressureDiscreteFunctionType PressureDiscreteFunctionType;
  #if PRESSURE_SPACE_TYPE == 2
  typedef typename Traits::PressureAdditionalDiscreteFunctionType PressureAdditionalDiscreteFunctionType;
  #endif
  typedef typename Traits::PressureDumpDiscreteFunctionType PressureDumpDiscreteFunctionType;
  typedef typename Traits::BulkDiscreteFunctionType BulkDiscreteFunctionType;
  typedef typename Traits::PhysicalCoefficientDiscreteFunctionType PhysicalCoefficientDiscreteFunctionType;
  typedef typename Traits::BulkDisplacementDiscreteFunctionType BulkDisplacementDiscreteFunctionType;
  typedef typename Traits::CurvatureDiscreteFunctionType CurvatureDiscreteFunctionType;
  typedef typename Traits::DisplacementDiscreteFunctionType DisplacementDiscreteFunctionType;
  typedef typename Traits::InterfaceDiscreteFunctionType InterfaceDiscreteFunctionType;
  // define tuple spaces
  typedef typename BulkDiscreteFunctionType::DiscreteFunctionSpaceType BulkDiscreteSpaceType;
  typedef typename InterfaceDiscreteFunctionType::DiscreteFunctionSpaceType InterfaceDiscreteSpaceType;
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
      bulkspace_=other.bulkspace_;
      #if PRESSURE_SPACE_TYPE == 2
      pressuredumpspace_=other.pressuredumpspace_;
      #endif
      bulkdisplacementspace_=other.bulkdisplacementspace_;
      #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
      rhospace_=other.rhospace_;
      #endif
      interfacespace_=other.interfacespace_;
      bulk_=other.bulk_;
      #if PRESSURE_SPACE_TYPE == 2
      pressuredump_=other.pressuredump_;
      #endif
      bulkdisplacement_=other.bulkdisplacement_;
      #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
      rho_=other.rho_;
      #endif
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
  const BulkGridType& bulkGrid() const
  {
    return meshmanager_.bulkGrid();
  }
  InterfaceGridType& interfaceGrid()
  {
    return meshmanager_.interfaceGrid();
  }
  const InterfaceGridType& interfaceGrid() const
  {
    return meshmanager_.interfaceGrid();
  }

  // get grid parts
  BulkGridPartType& bulkGridPart()
  {
    return meshmanager_.bulkGridPart();
  }
  const BulkGridPartType& bulkGridPart() const
  {
    return meshmanager_.bulkGridPart();
  }
  BulkInnerGridPartType& bulkInnerGridPart()
  {
    return meshmanager_.bulkInnerGridPart();
  }
  const BulkInnerGridPartType& bulkInnerGridPart() const
  {
    return meshmanager_.bulkInnerGridPart();
  }
  BulkOuterGridPartType& bulkOuterGridPart()
  {
    return meshmanager_.bulkOuterGridPart();
  }
  const BulkOuterGridPartType& bulkOuterGridPart() const
  {
    return meshmanager_.bulkOuterGridPart();
  }
  InterfaceGridPartType& interfaceGridPart()
  {
    return meshmanager_.interfaceGridPart();
  }
  const InterfaceGridPartType& interfaceGridPart() const
  {
    return meshmanager_.interfaceGridPart();
  }

  // get spaces
  const VelocityDiscreteSpaceType& velocitySpace() const
  {
    return bulkspace_->template subDiscreteFunctionSpace<0>();
  }
  const PressureDiscreteSpaceType& pressureSpace() const
  {
    return bulkspace_->template subDiscreteFunctionSpace<1>();
  }
  #if PRESSURE_SPACE_TYPE != 2
  const PressureDumpDiscreteSpaceType& pressureDumpSpace() const
  {
    return bulkspace_->template subDiscreteFunctionSpace<1>();
  }
  #else
  const PressureAdditionalDiscreteSpaceType& pressureAdditionalSpace() const
  {
    return bulkspace_->template subDiscreteFunctionSpace<2>();
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
  #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
  const PhysicalCoefficientDiscreteSpaceType& rhoSpace() const
  {
    return *rhospace_;
  }
  #endif
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
    return bulk_->template subDiscreteFunction<0>();
  }
  const VelocityDiscreteFunctionType& velocity() const
  {
    return bulk_->template subDiscreteFunction<0>();
  }
  PressureDiscreteFunctionType& pressure()
  {
    return bulk_->template subDiscreteFunction<1>();
  }
  const PressureDiscreteFunctionType& pressure() const
  {
    return bulk_->template subDiscreteFunction<1>();
  }
  #if PRESSURE_SPACE_TYPE != 2
  PressureDumpDiscreteFunctionType& pressureDump()
  {
    return bulk_->template subDiscreteFunction<1>();
  }
  const PressureDumpDiscreteFunctionType& pressureDump() const
  {
    return bulk_->template subDiscreteFunction<1>();
  }
  #else
  PressureAdditionalDiscreteFunctionType& pressureAdditional()
  {
    return bulk_->template subDiscreteFunction<2>();
  }
  const PressureAdditionalDiscreteFunctionType& pressureAdditional() const
  {
    return bulk_->template subDiscreteFunction<2>();
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
  #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
  PhysicalCoefficientDiscreteFunctionType& rho()
  {
    return *rho_;
  }
  const PhysicalCoefficientDiscreteFunctionType& rho() const
  {
    return *rho_;
  }
  #endif
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
    // reset pointers to avoid dangling references
    resetPointers();
    // create spaces
    bulkspace_=std::make_shared<BulkDiscreteSpaceType>(bulkGridPart());
    #if PRESSURE_SPACE_TYPE == 2
    pressuredumpspace_=std::make_shared<PressureDumpDiscreteSpaceType>(bulkGridPart());
    #endif
    bulkdisplacementspace_=std::make_shared<BulkDisplacementDiscreteSpaceType>(bulkGridPart());
    #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
    rhospace_=std::make_shared<PhysicalCoefficientDiscreteSpaceType>(bulkGridPart());
    #endif
    interfacespace_=std::make_shared<InterfaceDiscreteSpaceType>(interfaceGridPart());
    // create discrete functions
    bulk_=std::make_shared<BulkDiscreteFunctionType>("bulk solution",bulkSpace());
    velocity().name()="velocity";
    pressure().name()="pressure";
    #if PRESSURE_SPACE_TYPE == 2
    pressureAdditional().name()="additional pressure";
    pressuredump_=std::make_shared<PressureDumpDiscreteFunctionType>("pressure",pressureDumpSpace());
    #endif
    bulkdisplacement_=std::make_shared<BulkDisplacementDiscreteFunctionType>("bulk displacement",bulkDisplacementSpace());
    #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
    rho_=std::make_shared<PhysicalCoefficientDiscreteFunctionType>("rho",rhoSpace());
    #endif
    interface_=std::make_shared<InterfaceDiscreteFunctionType>("interface solution",interfaceSpace());
    curvature().name()="curvature";
    displacement().name()="displacement";
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
      bulkSolution().clear();
      bulkDisplacement().clear();
      #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
      rho().clear();
      #endif
      interfaceSolution().clear();
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
    Hybrid::forEach(typename BulkDiscreteFunctionType::Sequence{},
      [&](auto i){s<<" P"<<std::get<i>(*bulk_).space().order()<<" "<<std::get<i>(*bulk_).name()<<" -> "
        <<std::get<i>(*bulk_).size()<<" DOFs ";});
  }

  // print interface spaces info
  const void printInterfaceInfo(std::ostream& s=std::cout) const
  {
    Hybrid::forEach(typename BulkDiscreteFunctionType::Sequence{},
      [&](auto i){s<<" P"<<std::get<i>(*interface_).space().order()<<" "<<std::get<i>(*interface_).name()<<" -> "
        <<std::get<i>(*interface_).size()<<" DOFs ";});
  }

  // check pointers status
  void checkPointersStatus(const std::string& str="FluidState",std::ostream& s=std::cout) const
  {
    s<<std::endl;
    s<<str<<" number of pointers for each object (sequence = "<<sequence_<<") :"<<std::endl;
    s<<"BulkDiscreteSpace = "<<bulkspace_.use_count()<<std::endl;
    #if PRESSURE_SPACE_TYPE == 2
    s<<"PressureDumpDiscreteSpace = "<<pressuredumpspace_.use_count()<<std::endl;
    #endif
    s<<"BulkDisplacementDiscreteSpace = "<<bulkdisplacementspace_.use_count()<<std::endl;
    #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
    s<<"RhoDiscreteSpace = "<<rhospace_.use_count()<<std::endl;
    #endif
    s<<"InterfaceDiscreteSpace = "<<interfacespace_.use_count()<<std::endl;
    s<<"BulkDiscreteFunction = "<<bulk_.use_count()<<std::endl;
    #if PRESSURE_SPACE_TYPE == 2
    s<<"PressureDumpDiscreteFunction = "<<pressuredump_.use_count()<<std::endl;
    #endif
    s<<"BulkDisplacementDiscreteFunction = "<<bulkdisplacement_.use_count()<<std::endl;
    #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
    s<<"RhoDiscreteFunction = "<<rho_.use_count()<<std::endl;
    #endif
    s<<"InterfaceDiscreteFunction = "<<interface_.use_count()<<std::endl;
    s<<"BulkDataOutput = "<<bulkoutput_.use_count()<<std::endl;
    s<<"InterfaceDataOutput = "<<interfaceoutput_.use_count()<<std::endl;
    s<<std::endl;
  }

  private:
  CoupledMeshManagerType meshmanager_;
  std::shared_ptr<BulkDiscreteSpaceType> bulkspace_;
  #if PRESSURE_SPACE_TYPE == 2
  std::shared_ptr<PressureDumpDiscreteSpaceType> pressuredumpspace_;
  #endif
  std::shared_ptr<BulkDisplacementDiscreteSpaceType> bulkdisplacementspace_;
  #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
  std::shared_ptr<PhysicalCoefficientDiscreteSpaceType> rhospace_;
  #endif
  std::shared_ptr<InterfaceDiscreteSpaceType> interfacespace_;
  std::shared_ptr<BulkDiscreteFunctionType> bulk_;
  #if PRESSURE_SPACE_TYPE == 2
  std::shared_ptr<PressureDumpDiscreteFunctionType> pressuredump_;
  #endif
  std::shared_ptr<BulkDisplacementDiscreteFunctionType> bulkdisplacement_;
  #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
  std::shared_ptr<PhysicalCoefficientDiscreteFunctionType> rho_;
  #endif
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
    #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
    rho_.reset();
    #endif
    bulkdisplacement_.reset();
    #if PRESSURE_SPACE_TYPE == 2
    pressuredump_.reset();
    #endif
    bulk_.reset();
    interfacespace_.reset();
    #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
    rhospace_.reset();
    #endif
    bulkdisplacementspace_.reset();
    #if PRESSURE_SPACE_TYPE == 2
    pressuredumpspace_.reset();
    #endif
    bulkspace_.reset();
  }
};

}
}

#endif // DUNE_FEM_FLUIDSTATE_HH
