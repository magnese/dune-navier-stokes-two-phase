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
  std::cout<<std::endl<<"[DEBUG"<<info<<"]"<<std::endl;
  std::cout<<"\t Total mesh volume = "<<volumes[0]<<std::endl;
  std::cout<<"\t Entity min volume = "<<volumes[1]<<std::endl;
  std::cout<<"\t Entity max volume = "<<volumes[2]<<std::endl;
  std::cout<<"[DEBUG"<<info<<"]"<<std::endl<<std::endl;
  return volumes;
}

// dump interface volume
struct InterfaceVolumeInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  InterfaceVolumeInfo(unsigned int precision=6):
    BaseType("interface_volume",precision)
  {}

  using BaseType::add;

  template<typename CoupledMeshManagerType,typename TimeProviderType>
  void add(const CoupledMeshManagerType& meshManager,const TimeProviderType& timeProvider)
  {
    add(timeProvider.time(),meshManager.interfaceLength());
  }
};

// dump bulk inner volume
struct BulkInnerVolumeInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  BulkInnerVolumeInfo(unsigned int precision=6):
    BaseType("bulk_inner_volume",precision)
  {}

  using BaseType::add;

  template<typename CoupledMeshManagerType,typename TimeProviderType>
  void add(const CoupledMeshManagerType& meshManager,const TimeProviderType& timeProvider)
  {
    add(timeProvider.time(),meshManager.bulkInnerVolume());
  }
};

// dump normalized bulk inner volume
struct BulkNormalizedInnerVolumeInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  BulkNormalizedInnerVolumeInfo(unsigned int precision=6):
    BaseType("bulk_normalized_inner_volume",precision)
  {}

  using BaseType::add;
  using BaseType::isEmpty;
  using BaseType::firstValue;
  using BaseType::normalize;

  template<typename CoupledMeshManagerType,typename TimeProviderType>
  void add(const CoupledMeshManagerType& meshManager,const TimeProviderType& timeProvider)
  {
    add(timeProvider.time(),meshManager.bulkInnerVolume());
  }

  ~BulkNormalizedInnerVolumeInfo()
  {
    // normalize by intial volume
    if(!isEmpty())
      normalize(firstValue().second);
  }
};

// check if bulk and interface are consistent
template<typename InterfaceGridPartType,typename BulkGridPartType,typename BulkInterfaceGridMapperType>
bool isBulkConsistentWithInterface(const InterfaceGridPartType& interfaceGridPart,const BulkGridPartType& bulkGridPart,
                                   const BulkInterfaceGridMapperType& mapper)
{
  // perform an interface walkthrough
  const unsigned int bulkGriddim(BulkGridPartType::dimension);
  for(const auto& interfaceEntity:elements(interfaceGridPart))
  {
    // extract the corresponding bulk entity
    const auto interfaceIdx(interfaceGridPart.grid().leafIndexSet().index(interfaceEntity));
    const auto bulkEntity(bulkGridPart.grid().entity(mapper.entitySeedInterface2Bulk(interfaceIdx)));
    const auto& faceLocalIndex(mapper.faceLocalIdxInterface2Bulk(interfaceIdx));
    // create reference element
    const auto& refElement(ReferenceElements<typename BulkGridPartType::GridType::ctype,bulkGriddim>::general(bulkEntity.type()));
    // check inconsitency between bulk and interface
    for(auto i=decltype(bulkGriddim){0};i!=bulkGriddim;++i)
    {
      const auto bulkVtxLocalIndex(refElement.subEntity(faceLocalIndex,1,i,bulkGriddim));
      if(interfaceEntity.geometry().corner(i)!=bulkEntity.geometry().corner(bulkVtxLocalIndex))
      {
        std::cout<<std::endl<<"[DEBUG]"<<std::endl;
        std::cout<<"\t WARNING: Bulk and interface are NOT consistent!"<<std::endl;
        std::cout<<"[DEBUG]"<<std::endl<<std::endl;
        return false;
      }
    }
  }
  std::cout<<std::endl<<"[DEBUG]"<<std::endl;
  std::cout<<"\t Bulk and interface are consistent!"<<std::endl;
  std::cout<<"[DEBUG]"<<std::endl<<std::endl;
  return true;
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
  std::cout<<std::endl<<"[DEBUG]"<<std::endl;
  std::cout<<"\t |min("<<df.name()<<")| = "<<values[0]<<std::endl;
  std::cout<<"\t |max("<<df.name()<<")| = "<<values[1]<<std::endl;
  std::cout<<"\t |max("<<df.name()<<") - min("<<df.name()<<")| = "<<values[2]<<std::endl;
  std::cout<<"[DEBUG]"<<std::endl<<std::endl;
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
  std::cout<<std::endl<<"[DEBUG]"<<std::endl;
  std::cout<<"\t min("<<df.name()<<") = "<<values[0]<<std::endl;
  std::cout<<"\t max("<<df.name()<<") = "<<values[1]<<std::endl;
  std::cout<<"\t max("<<df.name()<<") - min("<<df.name()<<") = "<<values[2]<<std::endl;
  std::cout<<"[DEBUG]"<<std::endl<<std::endl;
  return values;
}

// dump function range
struct FunctionRangeInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  FunctionRangeInfo(unsigned int precision=6):
    BaseType("function_range",precision)
  {}

  using BaseType::add;

  template<typename DF,typename TimeProviderType>
  void add(const DF& df,const TimeProviderType& timeProvider)
  {
    std::array<double,2> values({std::numeric_limits<double>::max(),std::numeric_limits<double>::min()});
    for(const auto& dof:dofs(df))
    {
      values[0]=std::min(dof,values[0]);
      values[1]=std::max(dof,values[1]);
    }
    add(timeProvider.time(),values[1]-values[0]);
  }
};

// dump function max
struct FunctionMaxInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  FunctionMaxInfo(unsigned int precision=6):
    BaseType("function_max",precision)
  {}

  using BaseType::add;

  template<typename DF,typename TimeProviderType>
  void add(const DF& df,const TimeProviderType& timeProvider)
  {
    double value(std::numeric_limits<double>::min());
    for(const auto& dof:dofs(df))
      value=std::max(std::abs(dof),value);
    add(timeProvider.time(),value);
  }
};

// dump Tex log
template<typename TimerType,typename MeshManagerType>
void dumpTexLog(const std::vector<double>& errors,const TimerType& timer,const MeshManagerType& meshManager)
{
  std::ofstream file(Parameter::getValue<std::string>("fem.prefix",".")+"/log.tex");
  if(file.is_open())
  {
    file<<"\\documentclass[a4paper,11pt]{report}"<<std::endl<<std::endl;
    file<<"\\begin{document}"<<std::endl<<std::endl;
    file<<"$$\\mu_-="<<Parameter::getValue<double>("MuInner",1.0)<<"$$"<<std::endl;
    file<<"$$\\mu_+="<<Parameter::getValue<double>("MuOuter",1.0)<<"$$"<<std::endl;
    file<<"$$\\rho_-="<<Parameter::getValue<double>("RhoInner",1.0)<<"$$"<<std::endl;
    file<<"$$\\rho_+="<<Parameter::getValue<double>("RhoOuter",1.0)<<"$$"<<std::endl;
    file<<"$$\\gamma="<<Parameter::getValue<double>("Gamma",1.0)<<"$$"<<std::endl;
    file<<"$$C_s="<<Parameter::getValue<double>("CoeffSmoothing",1.0)<<"$$"<<std::endl;
    file<<"$$C_r="<<Parameter::getValue<double>("CoeffRemeshing",3.0)<<"$$"<<std::endl;
    file<<"$$\\tau="<<Parameter::getValue<double>("fem.timeprovider.fixedtimestep")<<"$$"<<std::endl;
    file<<"$$t\\in["<<Parameter::getValue<double>("fem.timeprovider.starttime",0)<<","
      <<Parameter::getValue<double>("EndTime",0.1)<<"]$$"<<std::endl<<std::endl;
    file<<"\\begin{table}"<<std::endl;
    file<<"\\center"<<std::endl;
    file<<"\\hspace*{-3.25cm}"<<std::endl;
    file<<"\\begin{tabular}{llllllllll}"<<std::endl;
    file<<"\\hline"<<std::endl;
    file<<"$K_\\Gamma$ & $\\|\\vec{X}-\\vec{x}\\|_{L^\\infty}$ & $\\|\\vec{U}-I^h_{2}\\vec{u}\\|_{L^2(\\Omega_T)}$"
      <<" & $\\|\\vec{U}-I^h_{2}\\vec{u}\\|_{H^1(\\Omega_T)}$ & $\\|\\vec{U}-I^h_{2}\\vec{u}\\|_{L^\\infty}$"
      <<" & $\\|P-p\\|_{L^2(\\Omega_T)}$ & $\\|P-I^h_{0}p\\|_{L^\\infty}$ & CPU[s] & $K_\\Omega^T$\\\\"<<std::endl;
    file<<"\\hline"<<std::endl;
    file<<std::setprecision(5);
    file<<std::scientific;
    file<<meshManager.interfaceGrid().size(0)<<" & ";
    const double nullTolerance(Parameter::getValue<double>("NullTolerance",1.e-12));
    for(const auto& err:errors)
      file<<(std::abs(err)<nullTolerance?0:err)<<" & ";
    file<<timer.elapsed()<<" & "<<meshManager.bulkGrid().size(0)<<"\\\\"<<std::endl;
    file<<"\\hline"<<std::endl;
    file<<"\\end{tabular}"<<std::endl;
    file<<"\\hspace*{-3.25cm}"<<std::endl;
    file<<"\\caption{Errors table.}"<<std::endl;
    file<<"\\end{table}"<<std::endl<<std::endl;
    file<<"\\end{document}"<<std::endl;
  }
}

// print discrete function boundary values with a certain boundaryID
template<typename DF,typename MeshManagerType>
void printDiscreteFunctionBoundaryValues(const DF& df,const MeshManagerType& meshManager,int boundaryID)
{
  std::cout<<std::endl<<"[DEBUG] "<<df.name()<<" boundary values with ID "<<boundaryID<<std::endl;
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
            std::cout<<std::endl;
          }
    }
  }
  std::cout<<"[DEBUG] "<<df.name()<<" boundary values with ID "<<boundaryID<<std::endl<<std::endl;
}

// dump triangles in gnuplot format
struct TrianglesDump:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  TrianglesDump(unsigned int precision=6):
    BaseType("triangles",precision)
  {}

  using BaseType::add;

  template<typename EntityType>
  void add(const EntityType& entity)
  {
    const auto& geo(entity.geometry());
    for(auto i=decltype(geo.corners()){0};i!=geo.corners();++i)
      add(geo.corner(i)[0],geo.corner(i)[1]);
    add(geo.corner(0)[0],geo.corner(0)[1],true);
  }
};


}
}

#endif // DUNE_FEM_MISCDEBUG_HH
