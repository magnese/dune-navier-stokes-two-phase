#ifndef DUNE_FEM_COUPLEDMESHMANAGER_HH
#define DUNE_FEM_COUPLEDMESHMANAGER_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/geometrygrid/grid.hh>
#include <dune/grid/io/file/gmshwriter.hh>
#include <dune/grid/utility/hostgridaccess.hh>
#include <dune/fem/gridpart/filter/basicfilterwrapper.hh>
#include <dune/fem/gridpart/filter/domainfilter.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/io/parameter.hh>

#include "vertexfunction.hh"
#include "gmshcompoundmanager.hh"

namespace Dune
{
namespace Fem
{

// volume criterion to trigger remeshing
class RemeshingVolumeCriterion
{
  public:
  typedef RemeshingVolumeCriterion ThisType;

  // constructor
  explicit RemeshingVolumeCriterion():
    coeff_(Parameter::getValidValue<double>("CoeffRemeshing",3.0,[](auto val){return val>=0.0;})),isenabled_(coeff_>0.0?true:false)
  {}

  RemeshingVolumeCriterion(const ThisType& )=default;

  ThisType& operator=(const ThisType& )
  {
    return *this;
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"Remeshing with volume criterion coefficient = "<<coeff_<<(isenabled_?"":" (WARNING: remesh disabled!)")<<"\n";
  }

  template<typename BulkGridPartType>
  bool remeshingIsNeeded(const BulkGridPartType& bulkGridPart) const
  {
    bool needed(false);
    if(isenabled_)
    {
      if(coeff_>1.0)
      {
        double maxVolume(std::numeric_limits<double>::min());
        double minVolume(std::numeric_limits<double>::max());
        for(const auto& entity:elements(bulkGridPart))
        {
          const auto volume(std::abs(entity.geometry().volume()));
          minVolume=std::min(volume,minVolume);
          maxVolume=std::max(volume,maxVolume);
          if(maxVolume/minVolume>coeff_)
          {
            std::cout<<"\nRemeshing needed.\n\n";
            needed=true;
            break;
          }
        }
      }
      else
      {
        std::cout<<"\nRemeshing performed at each time step.\n\n";
        needed=true;
      }
    }
    return needed;
  }

  private:
  const double coeff_;
  const bool isenabled_;
};

// angle criterion to trigger remeshing
class RemeshingAngleCriterion
{
  public:
  typedef RemeshingAngleCriterion ThisType;

  // constructor
  explicit RemeshingAngleCriterion():
    coeff_(Parameter::getValidValue<double>("CoeffRemeshing",20.0,[](auto val){return val>=0.0;})),isenabled_(coeff_>0.0?true:false)
  {}

  RemeshingAngleCriterion(const ThisType& )=default;

  ThisType& operator=(const ThisType& )
  {
    return *this;
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"Remeshing with angle criterion coefficient = "<<coeff_<<(isenabled_?"":" (WARNING: remesh disabled!)")<<"\n";
  }

  template<typename BulkGridPartType>
  bool remeshingIsNeeded(const BulkGridPartType& bulkGridPart) const
  {
    if(isenabled_)
    {
      constexpr auto rad2degCoeff(180.0/M_PI);
      typedef typename BulkGridPartType::GridType::template Codim<0>::Entity::Geometry::GlobalCoordinate GlobalCoordinateType;
      std::array<GlobalCoordinateType,BulkGridPartType::dimension+1> normals;
      for(const auto& entity:elements(bulkGridPart))
      {
        // store all the face normals
        unsigned int count(0);
        for(const auto& intersection:intersections(bulkGridPart,entity))
          normals[count++]=intersection.centerUnitOuterNormal();
        // compute all the possible angles between faces
        for(unsigned int i=0;i!=(normals.size()-1);++i)
          for(unsigned int j=(i+1);j!=normals.size();++j)
          {
            const auto angle(std::acos(-(normals[i]*normals[j]))*rad2degCoeff);
            if(angle<coeff_)
            {
              std::cout<<"\nRemeshing needed (found angle = "<<angle<<" ).\n\n";
              return true;
            }
          }
      }
    }
    return false;
  }

  private:
  const double coeff_;
  const bool isenabled_;
};

// coupled mesh manager
template<typename BulkHostGridImp,typename InterfaceHostGridImp,bool useCompoundManager,typename CharlengthPolicyImp=UniformCharlength,
         typename RemeshingCriterionImp=RemeshingVolumeCriterion>
class CoupledMeshManager
{
  public:
  // define host grids
  typedef BulkHostGridImp BulkHostGridType;
  typedef InterfaceHostGridImp InterfaceHostGridType;

  // define charlength policy and remeshing criterion
  typedef CharlengthPolicyImp CharlengthPolicyType;
  typedef RemeshingCriterionImp RemeshingCriterionType;

  typedef CoupledMeshManager<BulkHostGridType,InterfaceHostGridType,useCompoundManager,CharlengthPolicyType,RemeshingCriterionType>
    ThisType;

  // extract dimensions
  static constexpr unsigned int bulkGriddim=BulkHostGridType::dimension;
  static constexpr unsigned int interfaceGriddim=InterfaceHostGridType::dimension;
  static constexpr unsigned int worlddim=BulkHostGridType::dimensionworld;

  // grid factories types
  typedef GridFactory<BulkHostGridType> BulkHostGridFactoryType;
  typedef GridFactory<InterfaceHostGridType> InterfaceHostGridFactoryType;

  // define grids
  typedef GeometryGrid<BulkHostGridType,VertexFunction<BulkHostGridType>> BulkGridType;
  typedef GeometryGrid<InterfaceHostGridType,VertexFunction<InterfaceHostGridType>> InterfaceGridType;

  // define bulk and interface grid parts
  typedef LeafGridPart<BulkGridType> BulkGridPartType;
  typedef LeafGridPart<InterfaceGridType> InterfaceGridPartType;

  // define bulk indicator function
  typedef DomainFilter<BulkGridPartType,std::vector<int>> DomainFilterType;
  typedef BasicFilterWrapper<BulkGridPartType,DomainFilterType> IndicatorFunctionType;

  // define inner and outer grid parts
  #if PRESSURE_SPACE_TYPE == 0 || PRESSURE_SPACE_TYPE == 1 || PRESSURE_SPACE_TYPE == 2
  static constexpr bool useFilteredIndexSet=false;
  #elif PRESSURE_SPACE_TYPE == 3
  static constexpr bool useFilteredIndexSet=true;
  #endif
  typedef FilteredGridPart<BulkGridPartType,IndicatorFunctionType,useFilteredIndexSet> BulkInnerGridPartType;
  typedef FilteredGridPart<BulkGridPartType,IndicatorFunctionType,useFilteredIndexSet> BulkOuterGridPartType;

  // define interface entity to inner bulk entity mapper
  typedef typename BulkGridType::template Codim<0>::Entity::EntitySeed BulkEntitySeedType;
  typedef std::vector<std::pair<BulkEntitySeedType,int>> BulkInterfaceGridMapperType;
  typedef typename InterfaceGridType::template Codim<0>::Entity InterfaceEntityType;
  typedef typename BulkGridType::LeafIntersection BulkIntersectionType;

  // define bulk bounding box
  typedef std::pair<FieldVector<double,worlddim>,FieldVector<double,worlddim>> BulkBoundingBoxType;

  private:
  template<bool compound,std::size_t dim>
  struct GmshManagerSelector;

  template<std::size_t dim>
  struct GmshManagerSelector<true,dim>
  {
    typedef GMSHCompoundManager<dim,CharlengthPolicyType> GmshManagerType;
  };

  template<std::size_t dim>
  struct GmshManagerSelector<false,dim>
  {
    typedef GMSHSimpleManager GmshManagerType;
  };

  // define gmsh manager
  typedef typename GmshManagerSelector<useCompoundManager,worlddim>::GmshManagerType GmshManagerType;

  public:
  // constructor
  explicit CoupledMeshManager(int argc,char** argv,const GmshAlgorithmType& algorithm,bool verbosity,
                              bool checkEntityWithNoVerticesInDomain):
    manager_(argc,argv,algorithm,verbosity),sequence_(0),performentityverticescheck_(checkEntityWithNoVerticesInDomain)
  {
    init();
  }

  // copy constructor
  CoupledMeshManager(const ThisType& )=default;

  // copy assignament
  ThisType& operator=(const ThisType& other)
  {
    if(this!=&other)
    {
      resetPointers();
      boundaryids_=other.boundaryids_;
      elementsids_=other.elementsids_;
      manager_=other.manager_;
      bulkgrid_=other.bulkgrid_;
      bulkgridpart_=other.bulkgridpart_;
      bulkinnerindicator_=other.bulkinnerindicator_;
      bulkouterindicator_=other.bulkouterindicator_;
      bulkinnergridpart_=other.bulkinnergridpart_;
      bulkoutergridpart_=other.bulkoutergridpart_;
      interfacegrid_=other.interfacegrid_;
      interfacegridpart_=other.interfacegridpart_;
      mapper_=other.mapper_;
      remeshingcriterion_=other.remeshingcriterion_;
      sequence_=other.sequence_;
    }
    return *this;
  }

  void deepCopy(ThisType& other)
  {
    if(this!=&other)
    {
      std::cout<<"Deep copying of coupled mesh manager\n";
      // reset pointers to avoid dangling references
      resetPointers();
      // copy IDs
      boundaryids_=std::make_shared<std::vector<int>>(other.boundaryIDs());
      elementsids_=std::make_shared<std::vector<int>>(other.elementsIDs());
      // copy bulk grid
      std::vector<unsigned int> verticesMap(other.bulkGrid().size(bulkGriddim));
      unsigned int count(0);
      BulkHostGridFactoryType bulkHostGridFactory;
      for(const auto& vtx:vertices(other.bulkGridPart()))
      {
        bulkHostGridFactory.insertVertex(vtx.geometry().center());
        verticesMap[other.bulkGridPart().indexSet().index(vtx)]=(count++);
      }
      std::vector<unsigned int> entityConnectivity(bulkGriddim+1);
      for(const auto& entity:elements(other.bulkGridPart()))
      {
        for(auto vtxLocalIndex=decltype(entityConnectivity.size()){0};vtxLocalIndex!=entityConnectivity.size();++vtxLocalIndex)
          entityConnectivity[vtxLocalIndex]=verticesMap[other.bulkGridPart().indexSet().subIndex(entity,vtxLocalIndex,bulkGriddim)];
        bulkHostGridFactory.insertElement(entity.type(),entityConnectivity);
      }
      bulkgrid_=std::make_shared<BulkGridType>(bulkHostGridFactory.createGrid());
      printGridInfo(bulkGrid());
      // create bulk indicator functions and bulk grid parts
      bulkgridpart_=std::make_shared<BulkGridPartType>(bulkGrid());
      bulkinnerindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),bulkGridPart(),elementsIDs(),1);
      bulkouterindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),bulkGridPart(),elementsIDs(),2);
      bulkinnergridpart_=std::make_shared<BulkInnerGridPartType>(bulkGridPart(),bulkInnerIndicatorFunction());
      bulkoutergridpart_=std::make_shared<BulkOuterGridPartType>(bulkGridPart(),bulkOuterIndicatorFunction());
      // copy interface grid
      verticesMap.resize(other.interfaceGrid().size(interfaceGriddim));
      count=0;
      InterfaceHostGridFactoryType interfaceHostGridFactory;
      for(const auto& vtx:vertices(other.interfaceGridPart()))
      {
        interfaceHostGridFactory.insertVertex(vtx.geometry().center());
        verticesMap[other.interfaceGridPart().indexSet().index(vtx)]=(count++);
      }
      entityConnectivity.resize(interfaceGriddim+1);
      for(const auto& entity:elements(other.interfaceGridPart()))
      {
        for(auto vtxLocalIndex=decltype(entityConnectivity.size()){0};vtxLocalIndex!=entityConnectivity.size();++vtxLocalIndex)
          entityConnectivity[vtxLocalIndex]
            =verticesMap[other.interfaceGridPart().indexSet().subIndex(entity,vtxLocalIndex,interfaceGriddim)];
        interfaceHostGridFactory.insertElement(entity.type(),entityConnectivity);
      }
      interfacegrid_=std::make_shared<InterfaceGridType>(interfaceHostGridFactory.createGrid());
      printGridInfo(interfaceGrid());
      // create interface grid part
      interfacegridpart_=std::make_shared<InterfaceGridPartType>(interfaceGrid());
      // copy sequence, manager, remeshing criterion and mapper
      sequence_=other.sequence_;
      manager_=other.manager_;
      remeshingcriterion_=other.remeshingcriterion_;
      mapper_=std::make_shared<BulkInterfaceGridMapperType>(*(other.mapper_));
    }
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    manager_.printInfo(s);
    remeshingcriterion_.printInfo(s);
  }

  BulkGridType& bulkGrid()
  {
    return *bulkgrid_;
  }
  const BulkGridType& bulkGrid() const
  {
    return *bulkgrid_;
  }
  BulkGridPartType& bulkGridPart()
  {
    return *bulkgridpart_;
  }
  const BulkGridPartType& bulkGridPart() const
  {
    return *bulkgridpart_;
  }
  const IndicatorFunctionType& bulkInnerIndicatorFunction() const
  {
    return *bulkinnerindicator_;
  }
  const IndicatorFunctionType& bulkOuterIndicatorFunction() const
  {
    return *bulkouterindicator_;
  }
  BulkInnerGridPartType& bulkInnerGridPart()
  {
    return *bulkinnergridpart_;
  }
  const BulkInnerGridPartType& bulkInnerGridPart() const
  {
    return *bulkinnergridpart_;
  }
  BulkOuterGridPartType& bulkOuterGridPart()
  {
    return *bulkoutergridpart_;
  }
  const BulkOuterGridPartType& bulkOuterGridPart() const
  {
    return *bulkoutergridpart_;
  }
  InterfaceGridType& interfaceGrid()
  {
    return *interfacegrid_;
  }
  const InterfaceGridType& interfaceGrid() const
  {
    return *interfacegrid_;
  }
  InterfaceGridPartType& interfaceGridPart()
  {
    return *interfacegridpart_;
  }
  const InterfaceGridPartType& interfaceGridPart() const
  {
    return *interfacegridpart_;
  }
  GmshManagerType& manager()
  {
    return manager_;
  }
  const GmshManagerType& manager() const
  {
    return manager_;
  }

  // the sequence numebr increase each time the mesh changes
  unsigned int sequence() const
  {
    return sequence_;
  }

  // boundary ID of given intersection
  int intersectionID(const BulkIntersectionType& intersection) const
  {
    return boundaryIDs()[intersection.boundarySegmentIndex()];
  }

  // bulk inner entity intersection corresponding to the interface entity
  BulkIntersectionType correspondingInnerBulkIntersection(const InterfaceEntityType& interfaceEntity) const
  {
    const auto interfaceIndex(interfaceGrid().leafIndexSet().index(interfaceEntity));
    const auto bulkEntity(bulkGrid().entity((*mapper_)[interfaceIndex].first));
    const auto faceLocalIndex((*mapper_)[interfaceIndex].second);
    auto intersectionIt(bulkGridPart().ibegin(bulkEntity));
    while(intersectionIt->indexInInside()!=faceLocalIndex)
      ++intersectionIt;
    return *intersectionIt;
  }

  template<typename InterfaceFunctionType,typename BulkFunctionType>
  void setInterfaceDFInBulkDF(const InterfaceFunctionType& interfaceFunction,BulkFunctionType& bulkFunction) const
  {
    for(const auto& interfaceEntity:elements(interfaceGridPart()))
    {
      const auto localInterfaceFunction(interfaceFunction.localFunction(interfaceEntity));
      const auto intersection(correspondingInnerBulkIntersection(interfaceEntity));
      const auto bulkEntity(intersection.inside());
      auto localBulkFunction(bulkFunction.localFunction(bulkEntity));
      const auto faceLocalIndex(intersection.indexInInside());
      const auto numCorners(intersection.geometry().corners());
      const auto& refElement(ReferenceElements<typename BulkGridType::ctype,worlddim>::general(bulkEntity.type()));
      std::size_t offset(0);
      const auto bulkFirstVertex(
        HostGridAccess<BulkGridType>::hostEntity(bulkEntity).geometry().corner(refElement.subEntity(faceLocalIndex,1,0,worlddim)));
      while(HostGridAccess<InterfaceGridType>::hostEntity(interfaceEntity).geometry().corner(offset)!=bulkFirstVertex)
        ++offset;
      for(auto interfaceLocalIndex=decltype(numCorners){0};interfaceLocalIndex!=numCorners;++interfaceLocalIndex)
      {
        auto bulkRow(refElement.subEntity(faceLocalIndex,1,interfaceLocalIndex,worlddim)*worlddim);
        auto interfaceRow(((interfaceLocalIndex+offset)%worlddim)*worlddim);
        for(auto l=decltype(worlddim){0};l!=worlddim;++l,++interfaceRow,++bulkRow)
          localBulkFunction[bulkRow]=localInterfaceFunction[interfaceRow];
      }
    }
  }

  // list all boundary IDs
  std::vector<int> listUniqueIDs() const
  {
    std::unordered_set<int> uniqueIDs(10);
    for(const auto& boundaryID:boundaryIDs())
      uniqueIDs.insert(boundaryID);
    std::vector<int> vectorIDs(0);
    vectorIDs.reserve(uniqueIDs.size());
    for(const auto& boundaryID:uniqueIDs)
      vectorIDs.push_back(boundaryID);
    return vectorIDs;
  }

  bool remesh()
  {
    bool remeshPerformed(false);
    // check if the remeshing is supported
    if(GmshManagerType::remeshingSupported)
    {
      // check if the remeshing is necessary
      if(remeshingcriterion_.remeshingIsNeeded(bulkGridPart()))
      {
        // create timer
        Timer timer(false);
        timer.start();
        // reset pointers to avoid dangling references
        resetPointers(false);
        // create bulk grid and bulk grid part
        BulkHostGridFactoryType bulkHostGridFactory;
        boundaryids_=std::make_shared<std::vector<int>>();
        elementsids_=std::make_shared<std::vector<int>>();
        manager_.remesh(interfaceGrid(),bulkHostGridFactory,boundaryIDs(),elementsIDs());
        bulkgrid_=std::make_shared<BulkGridType>(bulkHostGridFactory.createGrid());
        printGridInfo(bulkGrid());
        bulkgridpart_=std::make_shared<BulkGridPartType>(bulkGrid());
        // reorder boundary IDs
        reorderBoundaryIDs(bulkHostGridFactory);
        // create bulk indicator functions and inner and outer grid parts
        bulkinnerindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),bulkGridPart(),elementsIDs(),1);
        bulkouterindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),bulkGridPart(),elementsIDs(),2);
        bulkinnergridpart_=std::make_shared<BulkInnerGridPartType>(bulkGridPart(),bulkInnerIndicatorFunction());
        bulkoutergridpart_=std::make_shared<BulkOuterGridPartType>(bulkGridPart(),bulkOuterIndicatorFunction());
        // create interface grid and interface grid part
        InterfaceHostGridFactoryType interfaceHostGridFactory;
        extractInterface(interfaceHostGridFactory);
        interfacegrid_=std::make_shared<InterfaceGridType>(interfaceHostGridFactory.createGrid());
        printGridInfo(interfaceGrid());
        interfacegridpart_=std::make_shared<InterfaceGridPartType>(interfaceGrid());
        // increase sequence number
        ++sequence_;
        // perform vertex check
        if(performentityverticescheck_)
          if(existEntityWithNoVerticesInDomain())
            DUNE_THROW(InvalidStateException,"ERROR: exists an entity with all the vertices on the boundary -> LBB not satisfied!");
        // print remesh time
        timer.stop();
        std::cout<<"Remeshing time: "<<timer.elapsed()<<" seconds.\n";
        remeshPerformed=true;
      }
    }
    return remeshPerformed;
  }

  double bulkVolume() const
  {
    double volume(0.0);
    for(const auto& entity:elements(bulkGridPart()))
      volume+=std::abs(entity.geometry().volume());
    return volume;
  }

  double bulkInnerVolume() const
  {
    double volume(0.0);
    for(const auto& entity:elements(bulkInnerGridPart()))
      volume+=std::abs(entity.geometry().volume());
    return volume;
  }

  double interfaceLength() const
  {
    double length(0.0);
    for(const auto& entity:elements(interfaceGridPart()))
      length+=std::abs(entity.geometry().volume());
    return length;
  }

  BulkBoundingBoxType bulkBoundingBox() const
  {
    FieldVector<double,worlddim> boundigBoxMin(std::numeric_limits<double>::max());
    FieldVector<double,worlddim> boundigBoxMax(std::numeric_limits<double>::min());
    for(const auto& vtx:vertices(bulkGridPart()))
    {
      const auto& pt(vtx.geometry().center());
      for(auto j=decltype(worlddim){0};j!=worlddim;++j)
      {
        boundigBoxMin[j]=std::min(boundigBoxMin[j],pt[j]);
        boundigBoxMax[j]=std::max(boundigBoxMax[j],pt[j]);
      }
    }
    return std::make_pair(std::move(boundigBoxMin),std::move(boundigBoxMax));
  }

  bool existEntityWithNoVerticesInDomain() const
  {
    for(const auto& entity:elements(bulkGridPart()))
      if(entity.hasBoundaryIntersections())
      {
        // count how many faces are boundary faces
        int count(0);
        for(const auto& intersection:intersections(bulkGridPart(),entity))
          if(intersection.boundary())
            ++count;
        // if the number of boundary faces is equal to worlddim it means that this entity has all vertices on the boundary
        if(count==worlddim)
          return true;
      }
    return false;
  }

  // check pointers status
  void checkPointersStatus(const std::string& str="CoupledMeshManager",std::ostream& s=std::cout) const
  {
    s<<"\n"<<str<<" number of pointers for each object (sequence = "<<sequence_<<") :\n";
    s<<"BoundaryIDs = "<<boundaryids_.use_count()<<"\n";
    s<<"ElementIDs = "<<elementsids_.use_count()<<"\n";
    s<<"BulkGrid = "<<bulkgrid_.use_count()<<"\n";
    s<<"BulkGridPart = "<<bulkgridpart_.use_count()<<"\n";
    s<<"BulkInnerIndicatorFunction = "<<bulkinnerindicator_.use_count()<<"\n";
    s<<"BulkOuterIndicatorFunction = "<<bulkouterindicator_.use_count()<<"\n";
    s<<"BulkInnerGridPart = "<<bulkinnergridpart_.use_count()<<"\n";
    s<<"BulkOuterGridPart = "<<bulkoutergridpart_.use_count()<<"\n";
    s<<"InterfaceGrid = "<<interfacegrid_.use_count()<<"\n";
    s<<"InterfaceGridPart = "<<interfacegridpart_.use_count()<<"\n";
    s<<"BulkInterfaceGridMapper = "<<mapper_.use_count()<<"\n\n";
  }

  // dump interface as msh file
  void dumpInterface(const std::string& fileName="interface.msh",int precision=16) const
  {
    GmshWriter<typename InterfaceGridType::LeafGridView> gmshWriter(interfaceGrid().leafGridView(),precision);
    gmshWriter.write(Parameter::getValue<std::string>("fem.prefix",".")+"/"+fileName,std::vector<int>(interfaceGrid().size(0),1));
  }

  // dump interface as msh file
  template<typename TimeProviderType>
  void dumpInterface(const TimeProviderType& timeProvider,int precision=16) const
  {
    dumpInterface("interface_"+std::to_string(timeProvider.time())+".msh",precision);
  }

  // dump bulk as msh file
  void dumpBulk(const std::string& fileName="bulk.msh",int precision=16) const
  {
    GmshWriter<typename BulkGridType::LeafGridView> gmshWriter(bulkGrid().leafGridView(),precision);
    gmshWriter.write(Parameter::getValue<std::string>("fem.prefix",".")+"/"+fileName,elementsIDs(),boundaryIDs());
  }

  // dump bulk as msh file
  template<typename TimeProviderType>
  void dumpBulk(const TimeProviderType& timeProvider,int precision=16) const
  {
    dumpBulk("bulk_"+std::to_string(timeProvider.time())+".msh",precision);
  }

  private:
  std::shared_ptr<std::vector<int>> boundaryids_;
  std::shared_ptr<std::vector<int>> elementsids_;
  GmshManagerType manager_;
  std::shared_ptr<BulkGridType> bulkgrid_;
  std::shared_ptr<BulkGridPartType> bulkgridpart_;
  std::shared_ptr<IndicatorFunctionType> bulkinnerindicator_;
  std::shared_ptr<IndicatorFunctionType> bulkouterindicator_;
  std::shared_ptr<BulkInnerGridPartType> bulkinnergridpart_;
  std::shared_ptr<BulkOuterGridPartType> bulkoutergridpart_;
  std::shared_ptr<InterfaceGridType> interfacegrid_;
  std::shared_ptr<InterfaceGridPartType> interfacegridpart_;
  std::shared_ptr<BulkInterfaceGridMapperType> mapper_;
  RemeshingCriterionType remeshingcriterion_;
  unsigned int sequence_;
  const bool performentityverticescheck_;

  std::vector<int>& boundaryIDs()
  {
    return *boundaryids_;
  }
  const std::vector<int>& boundaryIDs() const
  {
    return *boundaryids_;
  }
  std::vector<int>& elementsIDs()
  {
    return *elementsids_;
  }
  const std::vector<int>& elementsIDs() const
  {
    return *elementsids_;
  }

  // reset pointers to avoid dangling references
  void resetPointers(bool freeInterface=true)
  {
    mapper_.reset();
    interfacegridpart_.reset();
    if(freeInterface)
      interfacegrid_.reset();
    bulkoutergridpart_.reset();
    bulkinnergridpart_.reset();
    bulkouterindicator_.reset();
    bulkinnerindicator_.reset();
    bulkgridpart_.reset();
    bulkgrid_.reset();
    elementsids_.reset();
    boundaryids_.reset();
  }

  void init()
  {
    // reset pointers to avoid dangling references
    resetPointers();
    // create bulk grid and bulk grid part
    BulkHostGridFactoryType bulkHostGridFactory;
    boundaryids_=std::make_shared<std::vector<int>>();
    elementsids_=std::make_shared<std::vector<int>>();
    manager_.create(bulkHostGridFactory,boundaryIDs(),elementsIDs());
    bulkgrid_=std::make_shared<BulkGridType>(bulkHostGridFactory.createGrid());
    printGridInfo(bulkGrid());
    bulkgridpart_=std::make_shared<BulkGridPartType>(bulkGrid());
    // reorder boundary IDs
    reorderBoundaryIDs(bulkHostGridFactory);
    // create bulk indicator functions and inner and outer grid parts
    bulkinnerindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),bulkGridPart(),elementsIDs(),1);
    bulkouterindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),bulkGridPart(),elementsIDs(),2);
    bulkinnergridpart_=std::make_shared<BulkInnerGridPartType>(bulkGridPart(),bulkInnerIndicatorFunction());
    bulkoutergridpart_=std::make_shared<BulkOuterGridPartType>(bulkGridPart(),bulkOuterIndicatorFunction());
    // create interface grid and interface grid part
    InterfaceHostGridFactoryType interfaceHostGridFactory;
    extractInterface(interfaceHostGridFactory);
    interfacegrid_=std::make_shared<InterfaceGridType>(interfaceHostGridFactory.createGrid());
    printGridInfo(interfaceGrid());
    interfacegridpart_=std::make_shared<InterfaceGridPartType>(interfaceGrid());
    // increase sequence number
    ++sequence_;
    // perform vertex check
    if(performentityverticescheck_)
      if(existEntityWithNoVerticesInDomain())
        DUNE_THROW(InvalidStateException,"ERROR: exists an entity with all the vertices on the boundary -> LBB not satisfied!");
  }

  void extractInterface(InterfaceHostGridFactoryType& interfaceHostGridFactory)
  {
    // create new mapper
    mapper_=std::make_shared<BulkInterfaceGridMapperType>();
    std::size_t defaultValue(std::numeric_limits<std::size_t>::max());
    std::vector<std::size_t> vtxBulk2Interface(bulkGrid().size(worlddim),defaultValue);
    // create necessary variables
    std::size_t vtxInsertionCounter(0);
    GeometryType faceType(GeometryType::BasicType::simplex,interfaceGriddim);
    std::vector<unsigned int> faceConnectivity(bulkGriddim);
    // loop over bulk inner entities
    for(const auto& entity:elements(bulkInnerGridPart()))
      for(const auto& intersection:intersections(bulkInnerGridPart(),entity))
        if(intersection.neighbor())
          if(bulkOuterIndicatorFunction().contains(intersection.outside()))
          {
            const auto& refElement(ReferenceElements<typename BulkGridType::ctype,bulkGriddim>::general(entity.type()));
            const auto faceLocalIndex(intersection.indexInInside());
            for(auto i=decltype(bulkGriddim){0};i!=bulkGriddim;++i)
            {
              const auto vtxLocalIndex(refElement.subEntity(faceLocalIndex,1,i,bulkGriddim));
              const auto vtxGlobalIndex(bulkGrid().leafIndexSet().subIndex(entity,vtxLocalIndex,bulkGriddim));
              if(vtxBulk2Interface[vtxGlobalIndex]==defaultValue)
              {
                vtxBulk2Interface[vtxGlobalIndex]=vtxInsertionCounter++;
                interfaceHostGridFactory.insertVertex(entity.geometry().corner(vtxLocalIndex));
              }
              faceConnectivity[i]=vtxBulk2Interface[vtxGlobalIndex];
            }
            interfaceHostGridFactory.insertElement(faceType,faceConnectivity);
            mapper_->push_back(std::make_pair(entity.seed(),faceLocalIndex));
          }
  }

  void reorderBoundaryIDs(const BulkHostGridFactoryType& bulkHostGridFactory)
  {
    std::vector<int> tempIDs(boundaryIDs().size(),0);
    auto bulkHostLeafGridView(bulkGrid().hostGrid().leafGridView());
    for(const auto& entity:elements(bulkHostLeafGridView))
      for(const auto& intersection:intersections(bulkHostLeafGridView,entity))
        if(intersection.boundary())
          tempIDs[intersection.boundarySegmentIndex()]=boundaryIDs()[bulkHostGridFactory.insertionIndex(intersection)];
    boundaryIDs()=std::move(tempIDs);
  }

  template<typename GT>
  void printGridInfo(const GT& grid) const
  {
    std::cout<<"\nCreated grid (dimgrid = "<<GT::dimension<<") with "<<grid.size(0)<<" elements and "
      <<grid.size(GT::dimension)<<" vertices.\n\n";
  }
};

}
}

#endif // DUNE_FEM_COUPLEDMESHMANAGER_HH
