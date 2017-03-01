#ifndef DUNE_FEM_MISCDEBUG_HH
#define DUNE_FEM_MISCDEBUG_HH

// debug warning
#warning ("WARNING : you are using some debug functions!")

#include <iostream>
#include <array>
#include <algorithm>
#include <limits>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <utility>
#include <iomanip>

#include <dune/common/exceptions.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/rangegenerators.hh>

#include "gnuplotwriter.hh"

namespace Dune
{
namespace Fem
{

// compute total mesh volume, min and max entity volume
template<typename GridPartType>
std::array<double,3> meshVolumesInfo(const GridPartType& gridPart,const std::string& info="")
{
  std::array<double,3> volumes({0.0,std::numeric_limits<double>::max(),std::numeric_limits<double>::min()});
  for(const auto& entity:elements(gridPart))
  {
    const auto volume(std::abs(entity.geometry().volume()));
    volumes[0]+=volume;
    volumes[1]=std::min(volume,volumes[1]);
    volumes[2]=std::max(volume,volumes[2]);
  }
  std::cout<<"\n[DEBUG"<<info<<"]\n";
  std::cout<<"\t Total mesh volume = "<<volumes[0]<<"\n";
  std::cout<<"\t Entity min volume = "<<volumes[1]<<"\n";
  std::cout<<"\t Entity max volume = "<<volumes[2]<<"\n";
  std::cout<<"[DEBUG"<<info<<"]\n\n";
  return volumes;
}

// dump interface volume
struct InterfaceVolumeInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  InterfaceVolumeInfo(unsigned int precision=6):
    BaseType("interface_volume",precision)
  {}

  template<typename CoupledMeshManagerType,typename TimeProviderType>
  void add(const CoupledMeshManagerType& meshManager,const TimeProviderType& timeProvider)
  {
    BaseType::add(timeProvider.time(),meshManager.interfaceLength());
  }
};

// dump bulk inner volume
struct BulkInnerVolumeInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  BulkInnerVolumeInfo(unsigned int precision=6):
    BaseType("bulk_inner_volume",precision)
  {}

  template<typename CoupledMeshManagerType,typename TimeProviderType>
  void add(const CoupledMeshManagerType& meshManager,const TimeProviderType& timeProvider)
  {
    BaseType::add(timeProvider.time(),meshManager.bulkInnerVolume());
  }
};

// dump normalized bulk inner volume
struct BulkNormalizedInnerVolumeInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  BulkNormalizedInnerVolumeInfo(unsigned int precision=6):
    BaseType("bulk_normalized_inner_volume",precision)
  {}

  template<typename CoupledMeshManagerType,typename TimeProviderType>
  void add(const CoupledMeshManagerType& meshManager,const TimeProviderType& timeProvider)
  {
    BaseType::add(timeProvider.time(),meshManager.bulkInnerVolume());
  }

  ~BulkNormalizedInnerVolumeInfo()
  {
    // normalize by intial volume
    if(!BaseType::isEmpty())
      BaseType::normalize(BaseType::firstValue().second);
  }
};

// dump if remesh is triggered or not
struct RemeshInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  RemeshInfo(unsigned int precision=6):
    BaseType("remesh",precision)
  {}

  template<typename TimeProviderType>
  void add(bool remeshPerformed,const TimeProviderType& timeProvider)
  {
    if(remeshPerformed)
      BaseType::add(timeProvider.time(),1);
  }
};

// check if bulk and interface are consistent
template<typename FluidStateType>
void checkBulkInterfaceConsistency(FluidStateType& fluidState)
{
  // create bulk identity and interface identity
  fluidState.update();
  auto bulkIdentity(fluidState.bulkDisplacement());
  bulkIdentity.assign(fluidState.bulkGrid().coordFunction().discreteFunction());
  auto interfaceIdentity(fluidState.displacement());
  interfaceIdentity.assign(fluidState.interfaceGrid().coordFunction().discreteFunction());
  // create a copy of the bulk identity and set interface identity dofs in it
  auto bulkIdentityCopy(bulkIdentity);
  fluidState.meshManager().setInterfaceDFInBulkDF(interfaceIdentity,bulkIdentityCopy);
  // if they are different the setInterfaceDF2BulkDF is wrong
  bulkIdentityCopy-=bulkIdentity;
  for(const auto& dof:dofs(bulkIdentityCopy))
    if(dof!=0.0)
      DUNE_THROW(InvalidStateException,"ERROR: Bulk and interface are not consistent!");
}

// check if boundary normals are consistent
template<typename CoupledMeshManagerType>
void checkBoundaryNormalsConsistency(const CoupledMeshManagerType& meshManager,
  const typename CoupledMeshManagerType::BulkGridType::template Codim<0>::Entity::Geometry::GlobalCoordinate& center)
{
  bool isInitialized(false);
  double value(0);
  for(const auto& bulkEntity:elements(meshManager.bulkGridPart()))
    for(const auto& intersection:intersections(meshManager.bulkGridPart(),bulkEntity))
      if(intersection.boundary())
      {
        const auto normalVector(intersection.centerUnitOuterNormal());
        const auto positionVector(intersection.geometry().center()-center);
        const auto scalarProduct(normalVector.dot(positionVector));
        if(isInitialized)
        {
          if(value*scalarProduct<0)
            DUNE_THROW(InvalidStateException,"ERROR: Boundary normals are not consistent!");
        }
        else
        {
          value=scalarProduct;
          isInitialized=true;
        }
      }
}

// check if interface normals are consistent
template<typename CoupledMeshManagerType>
void checkInterfaceNormalsConsistency(const CoupledMeshManagerType& meshManager,
  const typename CoupledMeshManagerType::BulkGridType::template Codim<0>::Entity::Geometry::GlobalCoordinate& center)
{
  bool isInitialized(false);
  double value(0);
  for(const auto& interfaceEntity:elements(meshManager.interfaceGridPart()))
  {
    const auto intersection(meshManager.correspondingInnerBulkIntersection(interfaceEntity));
    const auto normalVector(intersection.centerUnitOuterNormal());
    const auto positionVector(interfaceEntity.geometry().center()-center);
    const auto scalarProduct(normalVector.dot(positionVector));
    if(isInitialized)
    {
      if(value*scalarProduct<0)
        DUNE_THROW(InvalidStateException,"ERROR: Interface normals are not consistent!");
    }
    else
    {
      value=scalarProduct;
      isInitialized=true;
    }
  }
}

// check abs function range
template<typename DF>
std::array<typename DF::RangeFieldType,3> checkFunctionAbsRange(const DF& df)
{
  typedef typename DF::RangeFieldType RangeFieldType;
  std::array<RangeFieldType,3> values({std::numeric_limits<RangeFieldType>::max(),std::numeric_limits<RangeFieldType>::min(),0.0});
  for(const auto& dof:dofs(df))
  {
    values[0]=std::min(std::abs(dof),values[0]);
    values[1]=std::max(std::abs(dof),values[1]);
  }
  values[2]=values[1]-values[0];
  std::cout<<"\n[DEBUG]\n";
  std::cout<<"\t |min("<<df.name()<<")| = "<<values[0]<<"\n";
  std::cout<<"\t |max("<<df.name()<<")| = "<<values[1]<<"\n";
  std::cout<<"\t |max("<<df.name()<<") - min("<<df.name()<<")| = "<<values[2]<<"\n";
  std::cout<<"[DEBUG]\n\n";
  return values;
}

// check function range
template<typename DF>
std::array<typename DF::RangeFieldType,3> checkFunctionRange(const DF& df)
{
  typedef typename DF::RangeFieldType RangeFieldType;
  std::array<RangeFieldType,3> values({std::numeric_limits<RangeFieldType>::max(),std::numeric_limits<RangeFieldType>::min(),0.0});
  for(const auto& dof:dofs(df))
  {
    values[0]=std::min(dof,values[0]);
    values[1]=std::max(dof,values[1]);
  }
  values[2]=values[1]-values[0];
  std::cout<<"\n[DEBUG]\n";
  std::cout<<"\t min("<<df.name()<<") = "<<values[0]<<"\n";
  std::cout<<"\t max("<<df.name()<<") = "<<values[1]<<"\n";
  std::cout<<"\t max("<<df.name()<<") - min("<<df.name()<<") = "<<values[2]<<"\n";
  std::cout<<"[DEBUG]\n\n";
  return values;
}

// dump function range
struct FunctionRangeInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  FunctionRangeInfo(unsigned int precision=6):
    BaseType("function_range",precision)
  {}

  template<typename DF,typename TimeProviderType>
  void add(const DF& df,const TimeProviderType& timeProvider)
  {
    std::array<double,2> values({std::numeric_limits<double>::max(),std::numeric_limits<double>::min()});
    for(const auto& dof:dofs(df))
    {
      values[0]=std::min(dof,values[0]);
      values[1]=std::max(dof,values[1]);
    }
    BaseType::add(timeProvider.time(),values[1]-values[0]);
  }
};

// dump function max
struct FunctionMaxInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  FunctionMaxInfo(unsigned int precision=6):
    BaseType("function_max",precision)
  {}

  template<typename DF,typename TimeProviderType>
  void add(const DF& df,const TimeProviderType& timeProvider)
  {
    double value(std::numeric_limits<double>::min());
    for(const auto& dof:dofs(df))
      value=std::max(std::abs(dof),value);
    BaseType::add(timeProvider.time(),value);
  }
};

// dump Tex log
template<typename TimerType,typename MeshManagerType>
void dumpTexLog(const std::vector<double>& errors,const TimerType& timer,const MeshManagerType& meshManager)
{
  std::ofstream file(Parameter::getValue<std::string>("fem.prefix",".")+"/log.tex");
  if(file.is_open())
  {
    file<<"\\documentclass[a4paper,11pt]{report}\n\n";
    file<<"\\begin{document}\n\n";
    file<<"$$\\mu_-="<<Parameter::getValue<double>("MuInner",1.0)<<"$$\n";
    file<<"$$\\mu_+="<<Parameter::getValue<double>("MuOuter",1.0)<<"$$\n";
    file<<"$$\\rho_-="<<Parameter::getValue<double>("RhoInner",1.0)<<"$$\n";
    file<<"$$\\rho_+="<<Parameter::getValue<double>("RhoOuter",1.0)<<"$$\n";
    file<<"$$\\gamma="<<Parameter::getValue<double>("Gamma",1.0)<<"$$\n";
    file<<"$$C_s="<<Parameter::getValue<double>("CoeffSmoothing",1.0)<<"$$\n";
    file<<"$$C_r="<<Parameter::getValue<double>("CoeffRemeshing",3.0)<<"$$\n";
    file<<"$$\\tau="<<Parameter::getValue<double>("fem.timeprovider.fixedtimestep")<<"$$\n";
    file<<"$$t\\in["<<Parameter::getValue<double>("fem.timeprovider.starttime",0)<<","
      <<Parameter::getValue<double>("EndTime",0.1)<<"]$$\n\n";
    file<<"\\begin{table}\n";
    file<<"\\center\n";
    file<<"\\hspace*{-3.25cm}\n";
    file<<"\\begin{tabular}{llllllllll}\n";
    file<<"\\hline\n";
    file<<"$K_\\Gamma$ & $\\|\\vec{X}-\\vec{x}\\|_{L^\\infty}$ & $\\|\\vec{U}-I^h_{2}\\vec{u}\\|_{L^2(\\Omega_T)}$"
      <<" & $\\|\\vec{U}-I^h_{2}\\vec{u}\\|_{H^1(\\Omega_T)}$ & $\\|\\vec{U}-I^h_{2}\\vec{u}\\|_{L^\\infty}$"
      <<" & $\\|P-p\\|_{L^2(\\Omega_T)}$ & $\\|P-I^h_{0}p\\|_{L^\\infty}$ & CPU[s] & $K_\\Omega^T$\\\\\n";
    file<<"\\hline\n";
    file<<std::setprecision(5);
    file<<std::scientific;
    file<<meshManager.interfaceGrid().size(0)<<" & ";
    const double nullTolerance(Parameter::getValue<double>("NullTolerance",1.e-12));
    for(const auto& err:errors)
      file<<(std::abs(err)<nullTolerance?0:err)<<" & ";
    file<<timer.elapsed()<<" & "<<meshManager.bulkGrid().size(0)<<"\\\\\n";
    file<<"\\hline\n";
    file<<"\\end{tabular}\n";
    file<<"\\hspace*{-3.25cm}\n";
    file<<"\\caption{Errors table.}\n";
    file<<"\\end{table}\n\n";
    file<<"\\end{document}\n";
  }
}

// print discrete function boundary values with a certain boundary ID
template<typename DF,typename MeshManagerType>
void printDiscreteFunctionBoundaryValues(const DF& df,const MeshManagerType& meshManager,int boundaryID)
{
  std::cout<<"\n[DEBUG] "<<df.name()<<" boundary values with ID "<<boundaryID<<"\n";
  const auto& gridPart(df.gridPart());
  for(const auto& entity:entities(df))
  {
    if(entity.hasBoundaryIntersections())
    {
      const auto& blockMapper(df.space().blockMapper());
      std::vector<std::size_t> globalIdxs(blockMapper.numDofs(entity));
      blockMapper.map(entity,globalIdxs);
      std::vector<bool> globalBlockDofsFilter(blockMapper.numDofs(entity));
      for(const auto& intersection:intersections(gridPart,entity))
        if(intersection.boundary())
          if(meshManager.intersectionID(intersection)==boundaryID)
          {
            blockMapper.onSubEntity(entity,intersection.indexInInside(),1,globalBlockDofsFilter);
            for(auto i=decltype(globalIdxs.size()){0};i!=globalIdxs.size();++i)
              if(globalBlockDofsFilter[i])
                std::cout<<"\t"<<*(df.block(globalIdxs[i]));
            std::cout<<"\n";
          }
    }
  }
  std::cout<<"[DEBUG] "<<df.name()<<" boundary values with ID "<<boundaryID<<"\n\n";
}

// dump triangles in gnuplot format
struct TrianglesDump:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  TrianglesDump(unsigned int precision=6):
    BaseType("triangles",precision)
  {}

  template<typename EntityType>
  void add(const EntityType& entity)
  {
    const auto& geo(entity.geometry());
    for(auto i=decltype(geo.corners()){0};i!=geo.corners();++i)
      BaseType::add(geo.corner(i)[0],geo.corner(i)[1]);
    BaseType::add(geo.corner(0)[0],geo.corner(0)[1],true);
  }
};


}
}

#endif // DUNE_FEM_MISCDEBUG_HH
