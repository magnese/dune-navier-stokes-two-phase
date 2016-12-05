#ifndef DUNE_FEM_COUPLEDMESHMANAGER_HH
#define DUNE_FEM_COUPLEDMESHMANAGER_HH

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/geometrygrid/grid.hh>
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

// mapper between bulk grid and interface grid
template<typename BulkGridImp,typename InterfaceGridImp>
class BulkInterfaceGridMapper
{
  public:
  typedef BulkGridImp BulkGridType;
  typedef InterfaceGridImp InterfaceGridType;
  typedef BulkInterfaceGridMapper<BulkGridType,InterfaceGridType> ThisType;
  typedef typename BulkGridType::template Codim<0>::Entity BulkEntityType;
  typedef typename BulkEntityType::EntitySeed BulkEntitySeedType;

  BulkInterfaceGridMapper():
    vtxbulk2interface_(0),vtxinterface2bulk_(0)
  {}

  BulkInterfaceGridMapper(const ThisType& )=default;
  ThisType& operator=(const ThisType& )=default;

  std::vector<std::size_t>& vtxBulk2Interface()
  {
    return vtxbulk2interface_;
  }

  const std::vector<std::size_t>& vtxBulk2Interface() const
  {
    return vtxbulk2interface_;
  }

  std::size_t& vtxBulk2Interface(std::size_t i)
  {
    return vtxbulk2interface_[i];
  }

  std::size_t vtxBulk2Interface(std::size_t i) const
  {
    return vtxbulk2interface_[i];
  }

  std::vector<std::size_t>& vtxInterface2Bulk()
  {
    return vtxinterface2bulk_;
  }

  const std::vector<std::size_t>& vtxInterface2Bulk() const
  {
    return vtxinterface2bulk_;
  }

  std::size_t& vtxInterface2Bulk(std::size_t i)
  {
    return vtxinterface2bulk_[i];
  }

  std::size_t& vtxInterface2Bulk(std::size_t i) const
  {
    return vtxinterface2bulk_[i];
  }

  void addBulkEntity2Mapping(const BulkEntityType& entity,unsigned int faceLocalIdx)
  {
     entityinterface2bulk_.push_back(std::make_pair(entity.seed(),faceLocalIdx));
  }

  const BulkEntitySeedType& entitySeedInterface2Bulk(std::size_t i) const
  {
    return entityinterface2bulk_[i].first;
  }

  unsigned int faceLocalIdxInterface2Bulk(std::size_t i) const
  {
    return entityinterface2bulk_[i].second;
  }

  void free()
  {
    vtxbulk2interface_.clear();
    vtxinterface2bulk_.clear();
    entityinterface2bulk_.clear();
  }

  template<typename InterfaceFunction,typename BulkFunction>
  void setInterfaceDFInBulkDF(const InterfaceFunction& interfaceFunction,BulkFunction& bulkFunction) const
  {
    constexpr std::size_t interfaceLocalBlockSize(InterfaceFunction::DiscreteFunctionSpaceType::localBlockSize);
    const auto interfaceNumBlocks(interfaceFunction.blocks());
    auto interfaceIt(interfaceFunction.dbegin());
    for(auto i=decltype(interfaceNumBlocks){0};i!=interfaceNumBlocks;++i)
    {
      const auto blockPos(vtxinterface2bulk_[i]);
      for(auto l=decltype(interfaceLocalBlockSize){0};l!=interfaceLocalBlockSize;++l,++interfaceIt)
        bulkFunction.dofVector()[blockPos][l]=*interfaceIt;
    }
  }

  private:
  std::vector<std::size_t> vtxbulk2interface_;
  std::vector<std::size_t> vtxinterface2bulk_;
  std::vector<std::pair<BulkEntitySeedType,unsigned int>> entityinterface2bulk_;
};

// volume criteria to trigger remeshing
class RemeshingVolumeCriteria
{
  public:
  typedef RemeshingVolumeCriteria ThisType;

  // constructor
  explicit RemeshingVolumeCriteria():
    coeff_(Parameter::getValue<double>("CoeffRemeshing",3.0))
  {}

  RemeshingVolumeCriteria(const ThisType& )=default;

  ThisType& operator=(const ThisType& )
  {
    return *this;
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"Remeshing Coefficient = "<<coeff_<<(coeff_<1.0?" (WARNING: remesh disabled!)":"")<<std::endl;
  }

  template<typename BulkGridPartType>
  bool remeshingIsNeeded(const BulkGridPartType& bulkGridPart) const
  {
    bool needed(false);
    if(coeff_>1.0)
    {
      double maxVolume(std::numeric_limits<double>::min());
      double minVolume(std::numeric_limits<double>::max());
      for(const auto& entity:elements(bulkGridPart))
      {
        const auto volume(std::abs(entity.geometry().volume()));
        minVolume=std::min(volume,minVolume);
        maxVolume=std::max(volume,maxVolume);
      }
      if(maxVolume/minVolume>coeff_)
      {
        std::cout<<std::endl<<"Remeshing needed (min bulk volume = "<<minVolume<<"; max bulk volume = "<<maxVolume<<")."
          <<std::endl<<std::endl;
        needed=true;
      }
    }
    if(coeff_==1.0)
    {
      std::cout<<std::endl<<"Remeshing performed at each time step."<<std::endl<<std::endl;
      needed=true;
    }
    return needed;
  }

  private:
  const double coeff_;
};

// coupled mesh manager
template<typename BulkHostGridImp,typename InterfaceHostGridImp,bool useCompoundManager,typename CharlengthPolicyImp=UniformCharlength,
         typename RemeshingCriteriaImp=RemeshingVolumeCriteria>
class CoupledMeshManager
{
  public:
  // define host grids
  typedef BulkHostGridImp BulkHostGridType;
  typedef InterfaceHostGridImp InterfaceHostGridType;

  // define charlength policy and remeshing criteria
  typedef CharlengthPolicyImp CharlengthPolicyType;
  typedef RemeshingCriteriaImp RemeshingCriteriaType;

  typedef CoupledMeshManager<BulkHostGridType,InterfaceHostGridType,useCompoundManager,CharlengthPolicyType,RemeshingCriteriaType> ThisType;

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
  typedef DomainFilter<BulkGridPartType,std::vector<int>> IndicatorFunctionType;

  // define inner and outer grid parts
  typedef FilteredGridPart<BulkGridPartType,IndicatorFunctionType,false> BulkInnerGridPartType;
  typedef FilteredGridPart<BulkGridPartType,IndicatorFunctionType,false> BulkOuterGridPartType;

  // define mapper
  typedef BulkInterfaceGridMapper<BulkGridType,InterfaceGridType> BulkInterfaceGridMapperType;

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
      remeshingcriteria_=other.remeshingcriteria_;
      sequence_=other.sequence_;
    }
    return *this;
  }

  void deepCopy(ThisType& other)
  {
    if(this!=&other)
    {
      std::cout<<"Deep copying of coupled mesh manager"<<std::endl;
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
      bulkinnerindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),elementsIDs(),1);
      bulkouterindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),elementsIDs(),2);
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
      // copy sequence, manager, remeshing criteria and mapper
      sequence_=other.sequence_;
      manager_=other.manager_;
      remeshingcriteria_=other.remeshingcriteria_;
      mapper_=std::make_shared<BulkInterfaceGridMapperType>(other.mapper());
    }
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    manager_.printInfo(s);
    remeshingcriteria_.printInfo(s);
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
  BulkInterfaceGridMapperType& mapper()
  {
    return *mapper_;
  }
  const BulkInterfaceGridMapperType& mapper() const
  {
    return *mapper_;
  }
  GmshManagerType& manager()
  {
    return manager_;
  }
  const GmshManagerType& manager() const
  {
    return manager_;
  }
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

  // the sequence numebr increase each time the mesh changes
  unsigned int sequence() const
  {
    return sequence_;
  }

  bool remesh()
  {
    bool remeshPerformed(false);
    // check if the remeshing is supported
    if(GmshManagerType::remeshingSupported)
    {
      // check if the remeshing is necessary
      if(remeshingcriteria_.remeshingIsNeeded(bulkGridPart()))
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
        bulkinnerindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),elementsIDs(),1);
        bulkouterindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),elementsIDs(),2);
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
        std::cout<<"Remeshing time: "<<timer.elapsed()<<" seconds."<<std::endl;
        remeshPerformed=true;
      }
    }
    return remeshPerformed;
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
    s<<std::endl;
    s<<str<<" number of pointers for each object (sequence = "<<sequence_<<") :"<<std::endl;
    s<<"BoundaryIDs = "<<boundaryids_.use_count()<<std::endl;
    s<<"ElementIDs = "<<elementsids_.use_count()<<std::endl;
    s<<"BulkGrid = "<<bulkgrid_.use_count()<<std::endl;
    s<<"BulkGridPart = "<<bulkgridpart_.use_count()<<std::endl;
    s<<"BulkInnerIndicatorFunction = "<<bulkinnerindicator_.use_count()<<std::endl;
    s<<"BulkOuterIndicatorFunction = "<<bulkouterindicator_.use_count()<<std::endl;
    s<<"BulkInnerGridPart = "<<bulkinnergridpart_.use_count()<<std::endl;
    s<<"BulkOuterGridPart = "<<bulkoutergridpart_.use_count()<<std::endl;
    s<<"InterfaceGrid = "<<interfacegrid_.use_count()<<std::endl;
    s<<"InterfaceGridPart = "<<interfacegridpart_.use_count()<<std::endl;
    s<<"BulkInterfaceGridMapper = "<<mapper_.use_count()<<std::endl;
    s<<std::endl;
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
  RemeshingCriteriaType remeshingcriteria_;
  unsigned int sequence_;
  const bool performentityverticescheck_;

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
    bulkinnerindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),elementsIDs(),1);
    bulkouterindicator_=std::make_shared<IndicatorFunctionType>(bulkGridPart(),elementsIDs(),2);
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
    mapper().vtxBulk2Interface().resize(bulkGrid().levelGridView(0).size(worlddim),defaultValue);
    // create necessary variables
    std::size_t vtxInsertionCounter(0);
    typename BulkGridType::template Codim<0>::Entity::Geometry::GlobalCoordinate vtx;
    std::array<std::size_t,bulkGriddim> vtxGlobalIndex;
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
              vtxGlobalIndex[i]=bulkGrid().leafIndexSet().subIndex(entity,vtxLocalIndex,bulkGriddim);
              if(mapper().vtxBulk2Interface(vtxGlobalIndex[i])==defaultValue)
              {
                mapper().vtxBulk2Interface(vtxGlobalIndex[i])=vtxInsertionCounter;
                vtx=entity.geometry().corner(vtxLocalIndex);
                interfaceHostGridFactory.insertVertex(vtx);
                mapper().vtxInterface2Bulk().push_back(vtxGlobalIndex[i]);
                ++vtxInsertionCounter;
              }
              faceConnectivity[i]=mapper().vtxBulk2Interface(vtxGlobalIndex[i]);
            }
            interfaceHostGridFactory.insertElement(faceType,faceConnectivity);
            mapper().addBulkEntity2Mapping(entity,faceLocalIndex);
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
    std::cout<<std::endl<<"Created grid (dimgrid = "<<GT::dimension<<") with "<<grid.size(0)<<" elements and "
      <<grid.size(GT::dimension)<<" vertices."<<std::endl<<std::endl;
  }
};

}
}

#endif // DUNE_FEM_COUPLEDMESHMANAGER_HH
