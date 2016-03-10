#ifndef DUNE_FEM_COUPLEDMESHMANAGER_HH
#define DUNE_FEM_COUPLEDMESHMANAGER_HH

// C++ includes
#include <iostream>
#include <limits>
#include <algorithm>
#include <vector>
#include <utility>
#include <memory>
#include <string>

// dune includes
#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/geometrygrid/grid.hh>
#include <dune/fem/io/parameter.hh>

// local includes
#include "indicatorfunction.hh"
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
  typedef typename BulkGridType::template Codim<0>::Entity BulkEntityType;
  typedef typename BulkEntityType::EntitySeed BulkEntitySeedType;

  BulkInterfaceGridMapper():
    vtxbulk2interface_(0),vtxinterface2bulk_(0)
  {}

  std::vector<std::size_t>& vtxBulk2Interface()
  {
    return vtxbulk2interface_;
  }

  const std::vector<std::size_t>& vtxBulk2Interface() const
  {
    return vtxbulk2interface_;
  }

  std::size_t& vtxBulk2Interface(const std::size_t& i)
  {
    return vtxbulk2interface_[i];
  }

  const std::size_t& vtxBulk2Interface(const std::size_t& i) const
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

  std::size_t& vtxInterface2Bulk(const std::size_t& i)
  {
    return vtxinterface2bulk_[i];
  }

  const std::size_t& vtxInterface2Bulk(const std::size_t& i) const
  {
    return vtxinterface2bulk_[i];
  }

  void addBulkEntity2Mapping(const BulkEntityType& entity,const unsigned int& faceLocalIdx)
  {
     entityinterface2bulk_.push_back(std::make_pair(entity.seed(),faceLocalIdx));
  }

  const BulkEntitySeedType& entitySeedInterface2Bulk(const std::size_t& i) const
  {
    return entityinterface2bulk_[i].first;
  }

  const unsigned int& faceLocalIdxInterface2Bulk(const std::size_t& i) const
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
  void addInterfaceDF2BulkDF(const InterfaceFunction& interfaceFunction,BulkFunction& bulkFunction) const
  {
    const auto interfaceLocalBlockSize(InterfaceFunction::DiscreteFunctionSpaceType::localBlockSize);
    const auto interfaceNumBlocks(interfaceFunction.blocks());
    auto interfaceIt(interfaceFunction.dbegin());
    for(auto i=0;i!=interfaceNumBlocks;++i)
    {
      const auto blockPos(vtxinterface2bulk_[i]);
      for(auto l=0;l!=interfaceLocalBlockSize;++l,++interfaceIt)
        bulkFunction.dofVector()[blockPos][l]+=*interfaceIt;
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
  // constructor
  explicit RemeshingVolumeCriteria():
    coeff_(Parameter::getValue<double>("CoeffRemeshing",3.0))
  {}

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"coeff_remesh = "<<coeff_<<(coeff_<1.0?" (WARNING: remesh disabled!)":"")<<std::endl;
  }

  template<typename GridType>
  bool remshingIsNeeded(const GridType& grid) const
  {
    bool needed(false);
    if(coeff_>1.0)
    {
      double maxVolume(std::numeric_limits<double>::min());
      double minVolume(std::numeric_limits<double>::max());
      auto leafGridView(grid.leafGridView());
      for(const auto& entity:elements(leafGridView))
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
      std::cout<<std::endl<<"Remeshing performed at each time step."<<std::endl<<std::endl;;
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
  static constexpr auto bulkGriddim=BulkHostGridType::dimension;
  static constexpr auto interfaceGriddim=InterfaceHostGridType::dimension;
  static constexpr auto worlddim=BulkHostGridType::dimensionworld;

  // grid factories types
  typedef GridFactory<BulkHostGridType> BulkHostGridFactoryType;
  typedef GridFactory<InterfaceHostGridType> InterfaceHostGridFactoryType;

  // define grids
  typedef GeometryGrid<BulkHostGridType,VertexFunction<BulkHostGridType>> BulkGridType;
  typedef GeometryGrid<InterfaceHostGridType,VertexFunction<InterfaceHostGridType>> InterfaceGridType;

  // define bulk indicator function
  typedef IndicatorFunction<BulkGridType> BulkIndicatorFunctionType;

  // define mapper
  typedef BulkInterfaceGridMapper<BulkGridType,InterfaceGridType> BulkInterfaceGridMapperType;

  //define bulk bounding box
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
  explicit CoupledMeshManager(int argc,char** argv,const GmshAlgorithmType& algorithm,const bool& verbosity,
                              const bool& checkEntityWithNoVerticesInDomain):
    manager_(argc,argv,algorithm,verbosity),sequence_(0),performentityverticescheck_(checkEntityWithNoVerticesInDomain)
  {
    init();
  }

  // copy constructor
  CoupledMeshManager(const ThisType& )=default;

  // copy assignament
  ThisType& operator=(const ThisType& )=default;

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
  const BulkIndicatorFunctionType& bulkIndicatorFunction() const
  {
    return *bulkindicator_;
  }
  InterfaceGridType& interfaceGrid()
  {
    return *interfacegrid_;
  }
  const InterfaceGridType& interfaceGrid() const
  {
    return *interfacegrid_;
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
  const unsigned int& sequence() const
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
      if(remeshingcriteria_.remshingIsNeeded(bulkGrid()))
      {
        // create timer
        Timer timer(false);
        timer.start();
        // create bulk grid
        BulkHostGridFactoryType bulkHostGridFactory;
        boundaryids_=std::make_shared<std::vector<int>>();
        elementsids_=std::make_shared<std::vector<int>>();
        manager_.remesh(interfaceGrid(),bulkHostGridFactory,boundaryIDs(),elementsIDs());
        bulkgrid_=std::make_shared<BulkGridType>(bulkHostGridFactory.createGrid());
        printGridInfo(bulkGrid());
        // reorder boundary IDs
        reorderBoundaryIDs(bulkHostGridFactory);
        // create bulk indicator function
        bulkindicator_=std::make_shared<BulkIndicatorFunctionType>(bulkGrid(),elementsIDs());
        // create interface grid
        InterfaceHostGridFactoryType interfaceHostGridFactory;
        extractInterface(interfaceHostGridFactory);
        interfacegrid_=std::make_shared<InterfaceGridType>(interfaceHostGridFactory.createGrid());
        printGridInfo(interfaceGrid());
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
    const auto leafGridView(bulkGrid().leafGridView());
    for(const auto& entity:elements(leafGridView))
      if(bulkIndicatorFunction().isInner(entity))
        volume+=std::abs(entity.geometry().volume());
    return volume;
  }

  double interfaceLength() const
  {
    double length(0.0);
    const auto leafGridView(interfaceGrid().leafGridView());
    for(const auto& entity:elements(leafGridView))
      length+=std::abs(entity.geometry().volume());
    return length;
  }

  BulkBoundingBoxType bulkBoundingBox() const
  {
    FieldVector<double,worlddim> boundigBoxMin(std::numeric_limits<double>::max());
    FieldVector<double,worlddim> boundigBoxMax(std::numeric_limits<double>::min());
    const auto& coords(bulkGrid().coordFunction().discreteFunction().dofVector());
    for(std::size_t i=0;i!=bulkGrid().size(0);++i)
      for(std::size_t j=0;j!=worlddim;++j)
      {
        boundigBoxMin[j]=std::min(boundigBoxMin[j],coords[i][j]);
        boundigBoxMax[j]=std::max(boundigBoxMax[j],coords[i][j]);
      }
    return std::make_pair(std::move(boundigBoxMin),std::move(boundigBoxMax));
  }

  bool existEntityWithNoVerticesInDomain() const
  {
    const auto gridLeafView(bulkGrid().leafGridView());
    for(const auto& entity:elements(gridLeafView))
      if(entity.hasBoundaryIntersections())
      {
        // count how many faces are boundary faces
        int count(0);
        for(const auto& intersection:intersections(gridLeafView,entity))
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
    s<<"BoundaryIDs = "<<boundaryids_.use_count()<<std::endl;;
    s<<"ElementIDs = "<<elementsids_.use_count()<<std::endl;
    s<<"BulkGrid = "<<bulkgrid_.use_count()<<std::endl;
    s<<"BulkIndicatorFunction = "<<bulkindicator_.use_count()<<std::endl;
    s<<"InterfaceGrid = "<<interfacegrid_.use_count()<<std::endl;
    s<<"BulkInterfaceGridMapper = "<<mapper_.use_count()<<std::endl;
    s<<std::endl;
  }

  private:
  std::shared_ptr<std::vector<int>> boundaryids_;
  std::shared_ptr<std::vector<int>> elementsids_;
  GmshManagerType manager_;
  std::shared_ptr<BulkGridType> bulkgrid_;
  std::shared_ptr<BulkIndicatorFunctionType> bulkindicator_;
  std::shared_ptr<InterfaceGridType> interfacegrid_;
  std::shared_ptr<BulkInterfaceGridMapperType> mapper_;
  RemeshingCriteriaType remeshingcriteria_;
  unsigned int sequence_;
  const bool performentityverticescheck_;

  void init()
  {
    // create bulk grid
    BulkHostGridFactoryType bulkHostGridFactory;
    boundaryids_=std::make_shared<std::vector<int>>();
    elementsids_=std::make_shared<std::vector<int>>();
    manager_.create(bulkHostGridFactory,boundaryIDs(),elementsIDs());
    bulkgrid_=std::make_shared<BulkGridType>(bulkHostGridFactory.createGrid());
    printGridInfo(bulkGrid());
    // reorder boundary IDs
    reorderBoundaryIDs(bulkHostGridFactory);
    // create bulk indicator function
    bulkindicator_=std::make_shared<BulkIndicatorFunctionType>(bulkGrid(),elementsIDs());
    // create interface grid
    InterfaceHostGridFactoryType interfaceHostGridFactory;
    extractInterface(interfaceHostGridFactory);
    interfacegrid_=std::make_shared<InterfaceGridType>(interfaceHostGridFactory.createGrid());
    printGridInfo(interfaceGrid());
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
    mapper().vtxBulk2Interface().resize(bulkGrid().levelGridView(0).size(worlddim),-1);
    // create necessary variables
    std::size_t vtxInsertionCounter(0);
    typename BulkGridType::template Codim<0>::Entity::Geometry::GlobalCoordinate vertex;
    std::array<std::size_t,bulkGriddim> vtxGlobalIndex;
    GeometryType faceType(GeometryType::BasicType::simplex,interfaceGriddim);
    std::vector<unsigned int> faceConnectivity(bulkGriddim);
    // loop over bulk entities
    auto bulkLeafGridView(bulkGrid().leafGridView());
    for(const auto& entity:elements(bulkLeafGridView))
    {
      const auto& refElement(ReferenceElements<typename BulkGridType::ctype,bulkGriddim>::general(entity.type()));
      if(bulkindicator_->isInner(entity))
      {
        for(const auto& intersection:intersections(bulkLeafGridView,entity))
        {
          if(intersection.neighbor())
          {
            const auto outsideIntersection(make_entity(intersection.outside()));
            if(!(bulkindicator_->isInner(outsideIntersection)))
            {
              const auto faceLocalIndex(intersection.indexInInside());
              for(auto i=0;i!=bulkGriddim;++i)
              {
                const auto vtxLocalIndex(refElement.subEntity(faceLocalIndex,1,i,bulkGriddim));
                vtxGlobalIndex[i]=bulkGrid().leafIndexSet().subIndex(entity,vtxLocalIndex,bulkGriddim);
                if(mapper().vtxBulk2Interface(vtxGlobalIndex[i])==-1)
                {
                  mapper().vtxBulk2Interface(vtxGlobalIndex[i])=vtxInsertionCounter;
                  vertex=entity.geometry().corner(vtxLocalIndex);
                  interfaceHostGridFactory.insertVertex(vertex);
                  mapper().vtxInterface2Bulk().push_back(vtxGlobalIndex[i]);
                  ++vtxInsertionCounter;
                }
                faceConnectivity[i]=mapper().vtxBulk2Interface(vtxGlobalIndex[i]);
              }
              interfaceHostGridFactory.insertElement(faceType,faceConnectivity);
              mapper().addBulkEntity2Mapping(entity,faceLocalIndex);
            }
          }
        }
      }
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
      <<grid.size(GT::dimension)<<" vertices."<<std::endl<<std::endl;;
  }
};

}
}

#endif // DUNE_FEM_COUPLEDMESHMANAGER_HH
