#ifndef DUNE_GMSHCOMPOUNDMANAGER_HH
#define DUNE_GMSHCOMPOUNDMANAGER_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <list>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>
#include <utility>

// HAVE_BLAS
#ifdef HAVE_BLAS
#define HAVE_BLAS_DUNE HAVE_BLAS
#undef HAVE_BLAS
#endif
// HAVE_LAPACK
#ifdef HAVE_LAPACK
#define HAVE_LAPACK_DUNE HAVE_LAPACK
#undef HAVE_LAPACK
#endif
// HAVE_METIS
#ifdef HAVE_METIS
#define HAVE_METIS_DUNE HAVE_METIS
#undef HAVE_METIS
#endif
// HAVE_GMP
#ifdef HAVE_GMP
#define HAVE_GMP_DUNE HAVE_GMP
#undef HAVE_GMP
#endif
// HAVE_PETSC
#ifdef HAVE_PETSC
#define HAVE_PETSC_DUNE HAVE_PETSC
#undef HAVE_PETSC
#endif

#include "Gmsh.h"
#include "GModel.h"
#include "MLine.h"
#include "MTriangle.h"
#include "MTetrahedron.h"

// HAVE_BLAS
#ifdef HAVE_BLAS
#undef HAVE_BLAS
#endif
#ifdef DUNE_HAVE_BLAS
#define HAVE_BLAS DUNE_HAVE_BLAS
#undef DUNE_HAVE_BLAS
#endif
// HAVE_LAPACK
#ifdef HAVE_LAPACK
#undef HAVE_LAPACK
#endif
#ifdef DUNE_HAVE_LAPACK
#define HAVE_LAPACK DUNE_HAVE_LAPACK
#undef DUNE_HAVE_LAPACK
#endif
// HAVE_METIS
#ifdef HAVE_METIS
#undef HAVE_METIS
#endif
#ifdef DUNE_HAVE_METIS
#define HAVE_METIS DUNE_HAVE_METIS
#undef DUNE_HAVE_METIS
#endif
// HAVE_GMP
#ifdef HAVE_GMP
#undef HAVE_GMP
#endif
#ifdef DUNE_HAVE_GMP
#define HAVE_GMP DUNE_HAVE_GMP
#undef DUNE_HAVE_GMP
#endif
// HAVE_PETSC
#ifdef HAVE_PETSC
#undef HAVE_PETSC
#endif
#ifdef DUNE_HAVE_PETSC
#define HAVE_PETSC DUNE_HAVE_PETSC
#undef DUNE_HAVE_PETSC
#endif

#include <dune/fem/io/io.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/gmshreader.hh>

namespace Dune
{

// fixed charlength policy
struct FixedCharlength
{
  typedef FixedCharlength ThisType;

  template<typename... Args>
  FixedCharlength(const Args&... )
  {}

  FixedCharlength(const ThisType& )=default;
  ThisType& operator=(const ThisType& )=default;

  static void printInfo(std::ostream& s=std::cout)
  {
    s<<"Remesh type = fixed charlength\n";
  }

  template<typename... Args>
  double operator()(const GVertex& vtx,const Args&...) const
  {
    return vtx.prescribedMeshSizeAtVertex();
  }
};

// uniform charlength policy
struct UniformCharlength
{
  typedef UniformCharlength ThisType;

  template<typename InterfaceGridType>
  UniformCharlength(const InterfaceGridType& interfaceGrid):
    charlength_(0)
  {
    // average volume
    for(const auto& entity:elements(interfaceGrid.leafGridView()))
      charlength_+=std::abs(entity.geometry().volume());
    charlength_/=static_cast<double>(interfaceGrid.size(0));
    // in 3D charlength is not the average volume but the edge length of an equilateral triangle with that area
    if(InterfaceGridType::dimensionworld==3)
      charlength_=2.0*std::pow(3.0,-0.25)*std::pow(charlength_,0.5);
  }

  UniformCharlength(double interfaceAverageVolume,unsigned int worlddim):
    charlength_(interfaceAverageVolume)
  {
    // in 3D charlength is not the average volume but the edge length of an equilateral triangle with that area
    if(worlddim==3)
      charlength_=2.0*std::pow(3.0,-0.25)*std::pow(charlength_,0.5);
  }

  UniformCharlength(const ThisType& )=default;
  ThisType& operator=(const ThisType& )=default;

  static void printInfo(std::ostream& s=std::cout)
  {
    s<<"Remesh type = uniform charlength\n";
  }

  template<typename... Args>
  double operator()(const Args&... ) const
  {
    return charlength_;
  }

  double charlength_;
};

// adaptive charlength policy
struct AdaptiveCharlength
{
  typedef AdaptiveCharlength ThisType;

  template<typename InterfaceGridType>
  AdaptiveCharlength(const InterfaceGridType& interfaceGrid):
    charlength_(0)
  {
    // average volume
    for(const auto& entity:elements(interfaceGrid.leafGridView()))
      charlength_+=std::abs(entity.geometry().volume());
    charlength_/=static_cast<double>(interfaceGrid.size(0));
    // in 3D charlength is not the average volume but the edge length of an equilateral triangle with that area
    if(InterfaceGridType::dimensionworld==3)
      charlength_=2.0*std::pow(3.0,-0.25)*std::pow(charlength_,0.5);
  }

  AdaptiveCharlength(double interfaceAverageVolume,unsigned int worlddim):
    charlength_(interfaceAverageVolume)
  {
    // in 3D charlength is not the average volume but the edge length of an equilateral triangle with that area
    if(worlddim==3)
      charlength_=2.0*std::pow(3.0,-0.25)*std::pow(charlength_,0.5);
  }

  AdaptiveCharlength(const ThisType& )=default;
  ThisType& operator=(const ThisType& )=default;

  static void printInfo(std::ostream& s=std::cout)
  {
    s<<"Remesh type = adaptive charlength\n";
  }

  template<typename VertexType,typename InterfaceGridType>
  double operator()(const VertexType& vtx,const InterfaceGridType& interfaceGrid) const
  {
    double dist(8.0*charlength_);
    for(const auto& vertex:vertices(interfaceGrid.leafGridView()))
    {
      FieldVector<double,3> temp{vtx.x(),vtx.y(),vtx.z()};
      for(std::size_t i=0;i!=InterfaceGridType::dimensionworld;i++)
        temp[i]-=vertex.geometry().center()[i];
      dist=std::min(dist,temp.two_norm());
    }
    return dist;
  }

  template<typename VertexType>
  double operator()(const VertexType& vtx,const std::shared_ptr<GModel>& interfaceModel) const
  {
    double dist(8.0*charlength_);
    for(auto vtxIt=interfaceModel->firstVertex();vtxIt!=interfaceModel->lastVertex();++vtxIt)
    {
      FieldVector<double,3> temp{vtx.x()-(*vtxIt)->x(),vtx.y()-(*vtxIt)->y(),vtx.z()-(*vtxIt)->z()};
      dist=std::min(dist,temp.two_norm());
    }
    return dist;
  }

  double charlength_;
};

// gmsh simple manager to use when the remeshing is not supported
class GMSHSimpleManager
{
  public:

  typedef GMSHSimpleManager ThisType;

  GMSHSimpleManager(int ,char** ):
    filename_(static_cast<std::string>(MSHFILESDIR)+Fem::Parameter::getValue<std::string>("CompoundFileName","compound.msh"))
  {}

  GMSHSimpleManager(const ThisType& )=default;

  ThisType& operator=(const ThisType& )
  {
    return *this;
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"Compound mesh filename: "<<filename_<<"\n";
    s<<"Remesh type = none\n";
  }

  template<typename BulkGridType,typename... Args>
  void create(GridFactory<BulkGridType>& bulkGridFactory,std::vector<int>& boundaryIDs,std::vector<int>& elementsIDs,
                     const Args&...)
  {
    GmshReader<BulkGridType>::read(bulkGridFactory,filename_,boundaryIDs,elementsIDs);
  }

  template<typename... Args>
  void remesh(const Args&...)
  {}

  static constexpr bool remeshingSupported=false;

  private:
  const std::string filename_;
};

// base class
template<unsigned int dim,typename CharlengthPolicyImp,typename Imp>
class GMSHCompoundManagerBase
{
  public:
  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"Domain mesh filename: "<<domainfilename_<<"\n";
    s<<"Interface mesh filename: "<<interfacefilename_<<"\n";
    if(hashole_)
      s<<"Hole mesh filename: "<<holefilename_<<"\n";
    CharlengthPolicyType::printInfo(s);
  }

  template<typename BulkGridType>
  void create(GridFactory<BulkGridType>& bulkGridFactory,std::vector<int>& boundaryIDs,std::vector<int>& elementsIDs)
  {
    // load domain
    domain()=std::make_shared<GModel>();
    domain()->setFactory("Gmsh");
    domain()->readGEO(domainfilename_);
    // load interface
    interface()=std::make_shared<GModel>();
    interface()->setFactory("Gmsh");
    bool readInterfaceFromGeo(interfacefilename_.find(".geo")!=std::string::npos);
    double interfaceAverageVolume(0);
    if(readInterfaceFromGeo)
      interface()->readGEO(interfacefilename_);
    else
    {
      interface()->readMSH(interfacefilename_);
      interfaceAverageVolume=imp().interfaceMesh2GModel(interface());
    }
    // load hole (if present)
    hole()=std::make_shared<GModel>();
    hole()->setFactory("Gmsh");
    if(hashole_)
      hole()->readGEO(holefilename_);
    if(readInterfaceFromGeo)
      imp().createCompoundGeo(FixedCharlength());
    else
    {
      CharlengthPolicyType charlengthPolicy(interfaceAverageVolume,worlddim);
      imp().createCompoundGeo([&](const GVertex& vtx){return charlengthPolicy(vtx,interface());});
    }
    compound()->mesh(worlddim);
    imp().bulkMesh2Dune(bulkGridFactory,boundaryIDs,elementsIDs);
  }

  static constexpr bool remeshingSupported=true;

  template<typename InterfaceGridType,typename BulkGridType>
  void remesh(const InterfaceGridType& interfaceGrid,GridFactory<BulkGridType>& bulkGridFactory,std::vector<int>& boundaryIDs,
              std::vector<int>& elementsIDs)
  {
    imp().interfaceDune2GModel(interfaceGrid);
    CharlengthPolicyType charlengthPolicy(interfaceGrid);
    imp().createCompoundGeo([&](const GVertex& vtx){return charlengthPolicy(vtx,interfaceGrid);});
    compound()->mesh(worlddim);
    imp().bulkMesh2Dune(bulkGridFactory,boundaryIDs,elementsIDs);
  }

  protected:
  typedef Imp Implementation;
  typedef CharlengthPolicyImp CharlengthPolicyType;
  typedef GMSHCompoundManagerBase<dim,CharlengthPolicyType,Implementation> ThisType;

  GMSHCompoundManagerBase(int argc,char** argv):
    domainfilename_(static_cast<std::string>(GEOFILESDIR)+Fem::Parameter::getValue<std::string>("DomainGeometry","domain.geo")),
    interfacefilename_(static_cast<std::string>(MSHFILESDIR)+Fem::Parameter::getValue<std::string>("InterfaceGeometry","interface.msh")),
    holefilename_(static_cast<std::string>(GEOFILESDIR)+Fem::Parameter::getValue<std::string>("HoleGeometry","")),gmodelptrs_(),
    hashole_(holefilename_!=static_cast<std::string>(GEOFILESDIR))
  {
    // init gmsh
    GmshSetOption("General","Terminal",1.);
    GmshSetOption("General","Verbosity",Fem::Parameter::getValue<double>("GmshVerbosity",0));
    GmshSetOption("Mesh","Algorithm",Fem::Parameter::getValue<double>("GmshAlgorithm",2));
  }

  GMSHCompoundManagerBase(const ThisType& )=default;

  ThisType& operator=(const ThisType& other)
  {
    gmodelptrs_=other.gmodelptrs_;
    hashole_=other.hashole_;
    return *this;
  }

  Implementation& imp()
  {
    return static_cast<Implementation&>(*this);
  }

  std::shared_ptr<GModel>& domain()
  {
    return gmodelptrs_[0];
  }
  std::shared_ptr<GModel>& interface()
  {
    return gmodelptrs_[1];
  }
  std::shared_ptr<GModel>& hole()
  {
    return gmodelptrs_[2];
  }
  std::shared_ptr<GModel>& compound()
  {
    return gmodelptrs_[3];
  }

  static constexpr unsigned int worlddim=dim;

  public:
  bool hasHole() const
  {
    return hashole_;
  }
  void writeInterfaceGeo(const std::string& fileName="interface.geo")
  {
    const std::string& path(Fem::Parameter::getValue<std::string>("fem.prefix","."));
    if(!Fem::directoryExists(path))
      Fem::createDirectory(path);
    interface()->writeGEO(path+"/"+fileName,true,false);
  }
  void writeInterfaceMsh(const std::string& fileName="interface.msh")
  {
    const std::string& path(Fem::Parameter::getValue<std::string>("fem.prefix","."));
    if(!Fem::directoryExists(path))
      Fem::createDirectory(path);
    interface()->writeMSH(path+"/"+fileName,2.2,false,false);
  }
  void writeCompoundGeo(const std::string& fileName="compound.geo")
  {
    const std::string& path(Fem::Parameter::getValue<std::string>("fem.prefix","."));
    if(!Fem::directoryExists(path))
      Fem::createDirectory(path);
    compound()->writeGEO(path+"/"+fileName,true,false);
  }
  void writeCompoundMsh(const std::string& fileName="compound.msh")
  {
    const std::string& path(Fem::Parameter::getValue<std::string>("fem.prefix","."));
    if(!Fem::directoryExists(path))
      Fem::createDirectory(path);
    compound()->writeMSH(path+"/"+fileName,2.2,false,false);
  }

  // check pointers status
  void checkPointersStatus(const std::string& str="GMSHCompoundManager",std::ostream& s=std::cout) const
  {
    s<<"\n"<<str<<" number of pointers for each object :\n";
    s<<"Domain = "<<domain().use_count()<<"\n";
    s<<"Interface = "<<interface().use_count()<<"\n";
    s<<"Hole = "<<hole().use_count()<<"\n";
    s<<"Compound = "<<compound().use_count()<<"\n\n";
  }

  private:
  const std::string domainfilename_;
  const std::string interfacefilename_;
  const std::string holefilename_;
  std::array<std::shared_ptr<GModel>,4> gmodelptrs_;
  bool hashole_;
};

// different specialization according to the dimension
template<unsigned int dim,typename CharlengthPolicyType=UniformCharlength>
class GMSHCompoundManager;

// specialization for worlddim = 2
template<typename CharlengthPolicyType>
class GMSHCompoundManager<2,CharlengthPolicyType>:
  public GMSHCompoundManagerBase<2,CharlengthPolicyType,GMSHCompoundManager<2,CharlengthPolicyType>>
{
  typedef GMSHCompoundManagerBase<2,CharlengthPolicyType,GMSHCompoundManager<2,CharlengthPolicyType>> BaseType;
  friend BaseType;
  using BaseType::compound;
  using BaseType::interface;
  using BaseType::domain;
  using BaseType::hole;
  using BaseType::hasHole;
  using BaseType::worlddim;

  public:
  template<typename... Args>
  GMSHCompoundManager(Args&&... args):
    BaseType(std::forward<Args>(args)...)
  {}

  private:
  template<typename CharlengthType>
  void createCompoundGeo(const CharlengthType& charlength)
  {
    compound()=std::make_shared<GModel>();
    compound()->setFactory("Gmsh");
    // add domain to compound gmodel
    std::vector<GEdge*> domainEdges(0);
    addGModelToCompound(domain(),domainEdges,charlength);
    // add interface to compound gmodel
    std::vector<GEdge*> interfaceEdges(0);
    addGModelToCompound(interface(),interfaceEdges,FixedCharlength());
    // add hole to compound gmodel (if present)
    std::vector<GEdge*> holeEdges(0);
    if(hasHole())
      addGModelToCompound(hole(),holeEdges,charlength);
    // add line loops and faces to compound gmodel
    std::vector<std::vector<GEdge*>> outerLineLoop{domainEdges,interfaceEdges};
    (compound()->addPlanarFace(outerLineLoop))->addPhysicalEntity(2);
    std::vector<std::vector<GEdge*>> innerLineLoop{interfaceEdges};
    if(hasHole())
      innerLineLoop.push_back(holeEdges);
    (compound()->addPlanarFace(innerLineLoop))->addPhysicalEntity(1);
  }

  template<typename CharlengthType>
  void addGModelToCompound(std::shared_ptr<GModel>& model,std::vector<GEdge*>& edges,const CharlengthType& charlength)
  {
    long int vtxCounter(0);
    std::vector<GVertex*> vertices(0);
    std::array<GVertex*,worlddim> vtxPtr;
    std::vector<long int> verticesMap(model->getNumVertices()+1,-1);
    for(auto it=model->firstEdge();it!=model->lastEdge();++it)
    {
      // get edge physical ID
      const auto physicalID(((*it)->getPhysicalEntities())[0]);
      // add first vertex
      vtxPtr[0]=(*it)->getBeginVertex();
      if(verticesMap[vtxPtr[0]->tag()]==-1)
      {
        vertices.push_back(compound()->addVertex(vtxPtr[0]->x(),vtxPtr[0]->y(),vtxPtr[0]->z(),charlength(*(vtxPtr[0]))));
        verticesMap[vtxPtr[0]->tag()]=vtxCounter;
        ++vtxCounter;
      }
      vtxPtr[0]=vertices[verticesMap[vtxPtr[0]->tag()]];
      // add last vertex
      vtxPtr[1]=(*it)->getEndVertex();
      if(verticesMap[vtxPtr[1]->tag()]==-1)
      {
        vertices.push_back(compound()->addVertex(vtxPtr[1]->x(),vtxPtr[1]->y(),vtxPtr[1]->z(),charlength(*(vtxPtr[1]))));
        verticesMap[vtxPtr[1]->tag()]=vtxCounter;
        ++vtxCounter;
      }
      vtxPtr[1]=vertices[verticesMap[vtxPtr[1]->tag()]];
      // add edge
      edges.push_back(compound()->addLine(vtxPtr[0],vtxPtr[1]));
      edges.back()->addPhysicalEntity(physicalID);
    }
  }

  double interfaceMesh2GModel(std::shared_ptr<GModel>& model)
  {
    // create new gmodel
    auto newGModel(std::make_shared<GModel>());
    newGModel->setFactory("Gmsh");
    // index all the mesh vertices in a continuous sequence starting at 1
    model->indexMeshVertices(true,0,true);
    long int vtxCounter(0);
    std::vector<GVertex*> vertices(0);
    std::vector<long int> verticesMap(model->getMaxVertexNumber()+1,0);
    constexpr double charlength(1000);
    // add vertices
    for(auto i=decltype(verticesMap.size()){1};i!=verticesMap.size();++i)
    {
      auto vtxPtr(model->getMeshVertexByTag(i));
      if(vtxPtr!=nullptr)
      {
        if(vtxPtr->getIndex()>(-1))
        {
          vertices.push_back(newGModel->addVertex(vtxPtr->x(),vtxPtr->y(),vtxPtr->z(),charlength));
          verticesMap[i]=vtxCounter;
          ++vtxCounter;
        }
      }
    }
    // add edges
    double lineTotalVolume(0);
    long int lineCounter(0);
    std::array<long int,worlddim> posVtx;
    constexpr int physicalID(1);
    for(auto edgeIt=model->firstEdge();edgeIt!=model->lastEdge();++edgeIt)
    {
      for(auto i=decltype((*edgeIt)->lines.size()){0};i!=(*edgeIt)->lines.size();++i)
      {
        auto linePtr((*edgeIt)->lines[i]);
        lineTotalVolume+=linePtr->getLength();
        ++lineCounter;
        posVtx[0]=verticesMap[linePtr->getVertex(0)->getNum()];
        posVtx[1]=verticesMap[linePtr->getVertex(1)->getNum()];
        (newGModel->addLine(vertices[posVtx[0]],vertices[posVtx[1]]))->addPhysicalEntity(physicalID);
      }
    }
    // assign new gmodel
    model=newGModel;
    return lineTotalVolume/static_cast<double>(lineCounter);
  }

  template<typename InterfaceGridType>
  void interfaceDune2GModel(const InterfaceGridType& interfaceGrid)
  {
    // create new gmodel
    interface()=std::make_shared<GModel>();
    interface()->setFactory("Gmsh");
    constexpr double charlength(1000);
    constexpr int physicalID(1);
    // insert all simplicies and all vertices of the interface in the GModel
    auto interfaceLeafGridView(interfaceGrid.leafGridView());
    std::vector<GVertex*> vertices(interfaceLeafGridView.size(worlddim-1),nullptr);
    std::array<long int,worlddim> posVtx;
    for(const auto& entity:elements(interfaceLeafGridView))
    {
      for(auto i=decltype(worlddim){0};i!=worlddim;++i)
      {
        posVtx[i]=interfaceGrid.leafIndexSet().subIndex(entity,i,worlddim-1);
        if(vertices[posVtx[i]]==nullptr)
        {
          const auto vtx(entity.geometry().corner(i));
          vertices[posVtx[i]]=interface()->addVertex(vtx[0],vtx[1],0.0,charlength);
        }
      }
      (interface()->addLine(vertices[posVtx[0]],vertices[posVtx[1]]))->addPhysicalEntity(physicalID);
    }
  }

  template<typename BulkGridType>
  void bulkMesh2Dune(GridFactory<BulkGridType>& gridFactory,std::vector<int>& boundaryIDs,std::vector<int>& elementsIDs)
  {
    // index all the mesh vertices in a continuous sequence starting at 1
    typename BulkGridType::template Codim<0>::Entity::Geometry::GlobalCoordinate vtx(0.0);
    compound()->indexMeshVertices(true,0,true);
    long int vtxCounter(0);
    std::vector<long int> verticesMap(compound()->getMaxVertexNumber()+1,0);
    // insert vertices in the grid factory
    for(auto i=decltype(verticesMap.size()){1};i!=verticesMap.size();++i)
    {
      auto vtxPtr(compound()->getMeshVertexByTag(i));
      if(vtxPtr!=nullptr)
      {
        if(vtxPtr->getIndex()>(-1))
        {
          vtx[0]=vtxPtr->x();
          vtx[1]=vtxPtr->y();
          gridFactory.insertVertex(vtx);
          verticesMap[i]=vtxCounter;
          ++vtxCounter;
        }
      }
    }
    // insert simplices in the grid factory
    GeometryType entityType(GeometryType::BasicType::simplex,worlddim);
    std::vector<unsigned int> entityConnectivity(worlddim+1);
    elementsIDs.clear();
    for(auto faceIt=compound()->firstFace();faceIt!=compound()->lastFace();++faceIt)
    {
      const auto entityID(((*faceIt)->getPhysicalEntities())[0]);
      for(auto i=decltype((*faceIt)->triangles.size()){0};i!=(*faceIt)->triangles.size();++i)
      {
        auto simplexPtr((*faceIt)->triangles[i]);
        for(auto j=decltype(worlddim){0};j!=(worlddim+1);++j)
          entityConnectivity[j]=verticesMap[simplexPtr->getVertex(j)->getNum()];
        gridFactory.insertElement(entityType,entityConnectivity);
        elementsIDs.push_back(entityID);
      }
    }
    // insert boundary faces in the grid factory
    std::vector<unsigned int> boundaryConnectivity(worlddim);
    boundaryIDs.clear();
    for(auto edgeIt=compound()->firstEdge();edgeIt!=compound()->lastEdge();++edgeIt)
    {
      const auto boundaryID(((*edgeIt)->getPhysicalEntities())[0]);
      if(boundaryID>1)
        for(auto i=decltype((*edgeIt)->lines.size()){0};i!=(*edgeIt)->lines.size();++i)
        {
          auto linePtr((*edgeIt)->lines[i]);
          for(auto j=decltype(worlddim){0};j!=worlddim;++j)
            boundaryConnectivity[j]=verticesMap[linePtr->getVertex(j)->getNum()];
          gridFactory.insertBoundarySegment(boundaryConnectivity);
          boundaryIDs.push_back(boundaryID);
        }
    }
  }
};

// specialization for worlddim = 3
template<typename CharlengthPolicyType>
class GMSHCompoundManager<3,CharlengthPolicyType>:
  public GMSHCompoundManagerBase<3,CharlengthPolicyType,GMSHCompoundManager<3,CharlengthPolicyType>>
{
  typedef GMSHCompoundManagerBase<3,CharlengthPolicyType,GMSHCompoundManager<3,CharlengthPolicyType>> BaseType;
  friend BaseType;
  using BaseType::compound;
  using BaseType::interface;
  using BaseType::domain;
  using BaseType::hole;
  using BaseType::hasHole;
  using BaseType::worlddim;

  public:
  template<typename... Args>
  GMSHCompoundManager(Args&&... args):
    BaseType(std::forward<Args>(args)...)
  {}

  private:
  template<typename CharlengthType>
  void createCompoundGeo(const CharlengthType& charlength)
  {
    compound()=std::make_shared<GModel>();
    compound()->setFactory("Gmsh");
    // add domain to compound gmodel
    std::vector<GFace*> domainFaces(0);
    addGModelToCompound(domain(),domainFaces,charlength);
    // add interface to compound gmodel
    std::vector<GFace*> interfaceFaces(0);
    addGModelToCompound(interface(),interfaceFaces,FixedCharlength());
    // add hole to compound gmodel (if present)
    std::vector<GFace*> holeFaces(0);
    if(hasHole())
      addGModelToCompound(hole(),holeFaces,charlength);
    // add surface loops and volumes to compound gmodel
    std::vector<std::vector<GFace*>> outerSurfaceLoop{domainFaces,interfaceFaces};
    (compound()->addVolume(outerSurfaceLoop))->addPhysicalEntity(2);
    std::vector<std::vector<GFace*>> innerSurfaceLoop{interfaceFaces};
    if(hasHole())
      innerSurfaceLoop.push_back(holeFaces);
    (compound()->addVolume(innerSurfaceLoop))->addPhysicalEntity(1);
  }

  template<typename CharlengthType>
  void addGModelToCompound(std::shared_ptr<GModel>& model,std::vector<GFace*>& faces,const CharlengthType& charlength)
  {
    std::vector<GVertex*> vertices(0);
    std::array<GVertex*,2> vtxPtr{nullptr,nullptr};
    std::vector<long int> verticesMap(model->getNumVertices()+1,-1);
    long int vtxCounter(0);
    std::map<long int,GEdge*> edgesMap;
    // loop over faces
    for(auto faceIt=model->firstFace();faceIt!=model->lastFace();++faceIt)
    {
      // get face physical ID
      const auto physicalID(((*faceIt)->getPhysicalEntities())[0]);
      auto edgesList((*faceIt)->edges());
      std::vector<GEdge*> edges(edgesList.size(),0);
      long int edgeCounter(0);
      // loop over edges
      for(auto& edge:edgesList)
      {
        auto edgeMapIt(edgesMap.find(edge->tag()));
        if(edgeMapIt==edgesMap.end())
        {
          // add first vertex
          vtxPtr[0]=edge->getBeginVertex();
          if(verticesMap[vtxPtr[0]->tag()]==-1)
          {
            vertices.push_back(compound()->addVertex(vtxPtr[0]->x(),vtxPtr[0]->y(),vtxPtr[0]->z(),charlength(*(vtxPtr[0]))));
            verticesMap[vtxPtr[0]->tag()]=vtxCounter;
            ++vtxCounter;
          }
          vtxPtr[0]=vertices[verticesMap[vtxPtr[0]->tag()]];
          // add last vertex
          vtxPtr[1]=edge->getEndVertex();
          if(verticesMap[vtxPtr[1]->tag()]==-1)
          {
            vertices.push_back(compound()->addVertex(vtxPtr[1]->x(),vtxPtr[1]->y(),vtxPtr[1]->z(),charlength(*(vtxPtr[1]))));
            verticesMap[vtxPtr[1]->tag()]=vtxCounter;
            ++vtxCounter;
          }
          vtxPtr[1]=vertices[verticesMap[vtxPtr[1]->tag()]];
          // add edge
          edges[edgeCounter]=compound()->addLine(vtxPtr[0],vtxPtr[1]);
          edgesMap.emplace(edge->tag(),edges[edgeCounter]);
        }
        else
          edges[edgeCounter]=edgeMapIt->second;
        ++edgeCounter;
      }
      // add lineloop
      std::vector<std::vector<GEdge*>> lineLoop{edges};
      faces.push_back(compound()->addPlanarFace(lineLoop));
      faces.back()->addPhysicalEntity(physicalID);
    }
  }

  double interfaceMesh2GModel(std::shared_ptr<GModel>& model)
  {
    // create new gmodel
    auto newGModel(std::make_shared<GModel>());
    newGModel->setFactory("Gmsh");
    // index all the mesh vertices in a continuous sequence starting at 1
    model->indexMeshVertices(true,0,true);
    long int vtxCounter(0);
    std::vector<GVertex*> vertices(0);
    typedef std::list<GEdge*> EdgeList;
    std::vector<std::pair<long int,EdgeList>> verticesMap(model->getMaxVertexNumber()+1,std::make_pair(0,EdgeList()));
    constexpr double charlength(1000);
    // add vertices
    for(auto i=decltype(verticesMap.size()){1};i!=verticesMap.size();++i)
    {
      auto vtxPtr(model->getMeshVertexByTag(i));
      if(vtxPtr!=nullptr)
      {
        if(vtxPtr->getIndex()>(-1))
        {
          vertices.push_back(newGModel->addVertex(vtxPtr->x(),vtxPtr->y(),vtxPtr->z(),charlength));
          verticesMap[i].first=vtxCounter;
          ++vtxCounter;
        }
      }
    }
    // add simplices
    double simplexTotalVolume(0);
    long int simplexCounter(0);
    std::array<long int,worlddim> idVtx;
    constexpr int physicalID(1);
    std::vector<GEdge*> simplexEdges(worlddim,nullptr);
    std::map<long int,GEdge*> edgesMap;
    for(auto faceIt=model->firstFace();faceIt!=model->lastFace();++faceIt)
    {
      for(auto i=decltype((*faceIt)->triangles.size()){0};i!=(*faceIt)->triangles.size();++i)
      {
        auto simplexPtr((*faceIt)->triangles[i]);
        simplexTotalVolume+=simplexPtr->getVolume();
        ++simplexCounter;
        // loop over edges
        for(auto l=decltype(worlddim){0};l!=worlddim;++l)
        {
          // get index of the 2 vertices
          for(int k=0;k!=2;++k)
            idVtx[k]=simplexPtr->getVertex((l+k)%worlddim)->getNum();
          // check if the edge has already been inserted
          auto edge0It(verticesMap[idVtx[0]].second.begin());
          const auto edge0ItEnd(verticesMap[idVtx[0]].second.end());
          const auto edge1ItEnd(verticesMap[idVtx[1]].second.end());
          bool found(false);
          while((edge0It!=edge0ItEnd)&&(!found))
          {
            for(auto edge1It=verticesMap[idVtx[1]].second.begin();(edge1It!=edge1ItEnd)&&(!found);++edge1It)
              if((*edge0It)->tag()==(*edge1It)->tag())
                found=true;
            if(!found)
              ++edge0It;
          }
          // use the already inserted edge or create the new one
          if(found)
            simplexEdges[l]=*edge0It;
          else
          {
            simplexEdges[l]=newGModel->addLine(vertices[verticesMap[idVtx[0]].first],vertices[verticesMap[idVtx[1]].first]);
            verticesMap[idVtx[0]].second.push_back(simplexEdges[l]);
            verticesMap[idVtx[1]].second.push_back(simplexEdges[l]);
          }
        }
        std::vector<std::vector<GEdge*>> edgeLoop{simplexEdges};
        (newGModel->addPlanarFace(edgeLoop))->addPhysicalEntity(physicalID);
      }
    }
    // assign new gmodel
    model=newGModel;
    return simplexTotalVolume/static_cast<double>(simplexCounter);
  }

  template<typename InterfaceGridType>
  void interfaceDune2GModel(const InterfaceGridType& interfaceGrid)
  {
    // create new gmodel
    interface()=std::make_shared<GModel>();
    interface()->setFactory("Gmsh");
    constexpr double charlength(1000);
    constexpr int physicalID(1);
    // insert all simplicies and all vertices of the interface in the GModel
    auto interfaceLeafGridView(interfaceGrid.leafGridView());
    std::vector<GVertex*> vertices(interfaceLeafGridView.size(worlddim-1),nullptr);
    std::vector<GEdge*> edges(interfaceLeafGridView.size(worlddim-2),nullptr);
    std::array<long int,worlddim> posVtx;
    std::array<long int,worlddim> posEdge;
    for(const auto& entity:elements(interfaceLeafGridView))
    {
      const auto refElement(referenceElement(entity.geometry()));
      for(auto i=decltype(worlddim){0};i!=worlddim;++i)
      {
        posVtx[i]=interfaceGrid.leafIndexSet().subIndex(entity,i,worlddim-1);
        if(vertices[posVtx[i]]==nullptr)
        {
          const auto vtx(entity.geometry().corner(i));
          vertices[posVtx[i]]=interface()->addVertex(vtx[0],vtx[1],vtx[2],charlength);
        }
      }
      for(auto i=decltype(worlddim){0};i!=worlddim;++i)
      {
        posEdge[i]=interfaceGrid.leafIndexSet().subIndex(entity,i,worlddim-2);
        if(edges[posEdge[i]]==nullptr)
          edges[posEdge[i]]=interface()->addLine(vertices[posVtx[refElement.subEntity(i,1,0,worlddim-1)]],
                                                 vertices[posVtx[refElement.subEntity(i,1,1,worlddim-1)]]);
      }
      std::vector<std::vector<GEdge*>> edgeLoop(1,{edges[posEdge[0]],edges[posEdge[1]],edges[posEdge[2]]});
      (interface()->addPlanarFace(edgeLoop))->addPhysicalEntity(physicalID);
    }
  }

  template<typename BulkGridType>
  void bulkMesh2Dune(GridFactory<BulkGridType>& gridFactory,std::vector<int>& boundaryIDs,std::vector<int>& elementsIDs)
  {
    // index all the mesh vertices in a continuous sequence starting at 1
    typename BulkGridType::template Codim<0>::Entity::Geometry::GlobalCoordinate vtx(0.0);
    compound()->indexMeshVertices(true,0,true);
    long int vtxCounter(0);
    std::vector<long int> verticesMap(compound()->getMaxVertexNumber()+1,0);
    // insert vertices in the grid factory
    for(auto i=decltype(verticesMap.size()){1};i!=verticesMap.size();++i)
    {
      auto vtxPtr(compound()->getMeshVertexByTag(i));
      if(vtxPtr!=nullptr)
      {
        if(vtxPtr->getIndex()>(-1))
        {
          vtx[0]=vtxPtr->x();
          vtx[1]=vtxPtr->y();
          vtx[2]=vtxPtr->z();
          gridFactory.insertVertex(vtx);
          verticesMap[i]=vtxCounter;
          ++vtxCounter;
        }
      }
    }
    // insert simplices in the grid factory
    GeometryType entityType(GeometryType::BasicType::simplex,worlddim);
    std::vector<unsigned int> entityConnectivity(worlddim+1);
    elementsIDs.clear();
    for(auto regionIt=compound()->firstRegion();regionIt!=compound()->lastRegion();++regionIt)
    {
      const auto entityID(((*regionIt)->getPhysicalEntities())[0]);
      for(auto i=decltype((*regionIt)->tetrahedra.size()){0};i!=(*regionIt)->tetrahedra.size();++i)
      {
        auto simplexPtr((*regionIt)->tetrahedra[i]);
        for(auto j=decltype(worlddim){0};j!=(worlddim+1);++j)
          entityConnectivity[j]=verticesMap[simplexPtr->getVertex(j)->getNum()];
        gridFactory.insertElement(entityType,entityConnectivity);
        elementsIDs.push_back(entityID);
      }
    }
    // insert boundary faces in the grid factory
    std::vector<unsigned int> boundaryConnectivity(worlddim);
    boundaryIDs.clear();
    for(auto faceIt=compound()->firstFace();faceIt!=compound()->lastFace();++faceIt)
    {
      const auto boundaryID(((*faceIt)->getPhysicalEntities())[0]);
      if(boundaryID>1)
        for(auto i=decltype((*faceIt)->triangles.size()){0};i!=(*faceIt)->triangles.size();++i)
        {
          auto trianglePtr((*faceIt)->triangles[i]);
          for(auto j=decltype(worlddim){0};j!=worlddim;++j)
            boundaryConnectivity[j]=verticesMap[trianglePtr->getVertex(j)->getNum()];
          gridFactory.insertBoundarySegment(boundaryConnectivity);
          boundaryIDs.push_back(boundaryID);
        }
    }
  }
};

}
#endif // DUNE_GMSHCOMPOUNDMANAGER_HH
