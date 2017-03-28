#ifndef DUNE_FEM_BOUNDARYCONDITION_HH
#define DUNE_FEM_BOUNDARYCONDITION_HH

#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>

#include "extendedtuplediscretefunction.hh"
#include "gridparthelper.hh"

namespace Dune
{
namespace Fem
{

// struct to store block IDs when overlapping IDs can NOT be handled using only one of the IDs involved
class ListBlockID
{
  public:
  typedef std::list<int> value_type;

  ListBlockID(int val=-1):
    values_(1,val)
  {}
  void clear()
  {
    values_.clear();
  }
  void operator=(int value)
  {
    if(values_.front()==-1)
      values_.front()=value;
    else
      values_.push_back(value);
  }
  const value_type& get() const
  {
    return values_;
  }
  private:
  value_type values_;
};

// BC  type
enum BCEnumType{dirichlet,neumann,robin,freeslip};

// BC interface
template<typename DiscreteSpaceImp,typename CoupledMeshManagerImp,typename BlockIDImp,typename BCImp>
class BoundaryCondition
{
  public:
  typedef DiscreteSpaceImp DiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef BlockIDImp BlockIDType;
  typedef BCImp BCImplementation;
  typedef BoundaryCondition<DiscreteSpaceType,CoupledMeshManagerType,BlockIDType,BCImplementation> ThisType;

  typedef typename DiscreteSpaceType::DomainType DomainType;
  typedef typename DiscreteSpaceType::RangeType RangeType;
  typedef typename DiscreteSpaceType::EntityType EntityType;
  typedef typename DiscreteSpaceType::GridPartType::IntersectionIteratorType::Intersection IntersectionType;

  typedef LocalAnalyticalFunctionBinder<DiscreteSpaceType> LocalAnalyticalFunctionType;
  typedef std::map<int,LocalAnalyticalFunctionType> FunctionMapType;
  typedef LocalFunctionAdapter<LocalAnalyticalFunctionType> AdaptedDiscreteFunctionType;
  typedef std::map<int,AdaptedDiscreteFunctionType> AdaptedFunctionMapType;

  typedef DynamicVector<typename AdaptedDiscreteFunctionType::RangeFieldType> LocalBoundaryDOFsType;

  const DiscreteSpaceType& space() const
  {
    return *space_;
  }

  template<typename TimeProviderType,typename RHSType,typename... OperatorsType>
  void apply(const TimeProviderType& timeProvider,RHSType& rhs,OperatorsType&... operators)
  {
    update();
    updateDOFs();
    for(const auto& entity:space())
      asImp().setDOFs(entity,timeProvider,rhs,operators...);
  }

  RangeType evaluateBoundaryFunction(const DomainType& x,double t,const EntityType& entity,int boundaryID) const
  {
    auto gIt(g_.find(boundaryID));
    if(gIt==g_.end())
      DUNE_THROW(RangeError,"ERROR: boundary ID not found in BC!");
    auto& g(gIt->second);
    g.init(entity);
    g.initialize(t,t);
    RangeType ret;
    g.evaluate(x,ret);
    return ret;
  }

  RangeType evaluateBoundaryFunction(const DomainType& x,double t,const IntersectionType& intersection) const
  {
    return evaluateBoundaryFunction(x,t,intersection.inside(),meshmanager_.intersectionID(intersection));
  }

  LocalBoundaryDOFsType localBoundaryDOFs(double t,const EntityType& entity,int boundaryID) const
  {
    auto gAdaptedIt(gadapted_.find(boundaryID));
    if(gAdaptedIt==gadapted_.end())
      DUNE_THROW(RangeError,"ERROR: boundary ID not found in BC!");
    auto& gAdapted(gAdaptedIt->second);
    gAdapted.initialize(t,t);
    const auto interpolation(space().interpolation(entity));
    LocalBoundaryDOFsType localDOFs(space().basisFunctionSet(entity).size());
    interpolation(gAdapted.localFunction(entity),localDOFs);
    return localDOFs;
  }

  LocalBoundaryDOFsType localBoundaryDOFs(double t,const IntersectionType& intersection) const
  {
    return localBoundaryDOFs(t,intersection.inside(),meshmanager_.intersectionID(intersection));
  }

  template<typename FilteredIntersectionType>
  LocalBoundaryDOFsType localBoundaryDOFs(double t,const FilteredIntersectionType& intersection) const
  {
    return localBoundaryDOFs(t,intersection.inside(),meshmanager_.intersectionID(intersection.hostIntersection()));
  }

  // update space and gadapt
  void update()
  {
    // check if the mesh is changed
    if(sequence_!=meshmanager_.sequence())
    {
      // create space
      space_=std::make_unique<DiscreteSpaceType>(meshmanager_.bulkGridPart());
      // create adapted discrete boundary function
      gadapted_.clear();
      for(auto& g:g_)
        gadapted_.emplace(g.first,AdaptedDiscreteFunctionType("adapted function",g.second,meshmanager_.bulkGridPart(),
          space_->order()));
      // update sequence number
      sequence_=meshmanager_.sequence();
    }
  }

  BoundaryCondition()=delete;
  BoundaryCondition(const ThisType&)=delete;
  ThisType& operator=(const ThisType&)=delete;
  BoundaryCondition(ThisType&&)=default;
  ThisType& operator=(ThisType&&)=default;

  protected:
  CoupledMeshManagerType& meshmanager_;
  std::unique_ptr<DiscreteSpaceType> space_;
  const BCEnumType bctype_;
  mutable FunctionMapType g_;
  mutable AdaptedFunctionMapType gadapted_;
  std::vector<BlockIDType> blocksIDs_;
  unsigned int sequence_;
  unsigned int sequencedof_;

  BoundaryCondition(CoupledMeshManagerType& meshManager,BCEnumType bctype):
    meshmanager_(meshManager),bctype_(bctype),g_(),gadapted_(),blocksIDs_(0),sequence_(0),sequencedof_(0)
  {}

  BCImplementation& asImp()
  {
    return static_cast<BCImplementation&>(*this);
  }

  const BCImplementation& asImp() const
  {
    return static_cast<const BCImplementation&>(*this);
  }

  // update the blocks IDs when the mesh changes
  void updateDOFs()
  {
    // check if the mesh is changed
    if(sequencedof_!=meshmanager_.sequence())
    {
      // set blocksIDs
      blocksIDs_.clear();
      blocksIDs_.resize(space().blockMapper().size(),-1);
      for(const auto& entity:space())
      {
        if(entity.hasBoundaryIntersections())
          setBlocksIDs(entity);
      }
      // update sequence number
      sequencedof_=meshmanager_.sequence();
    }
  }

  // set the correct IDs for each block
  void setBlocksIDs(const EntityType& entity)
  {
    const auto& blockMapper(space().blockMapper());
    std::vector<std::size_t> globalIdxs(blockMapper.numDofs(entity));
    blockMapper.map(entity,globalIdxs);
    std::vector<bool> globalBlockDofsFilter(blockMapper.numDofs(entity));

    for(const auto& intersection:intersections(meshmanager_.bulkGridPart(),entity))
    {
      if(intersection.boundary())
      {
        const auto& boundaryID(meshmanager_.intersectionID(intersection));
        if(g_.find(boundaryID)!=g_.end())
        {
          const auto faceLocalIdx(intersection.indexInInside());
          blockMapper.onSubEntity(entity,faceLocalIdx,1,globalBlockDofsFilter);
          for(auto i=decltype(globalIdxs.size()){0};i!=globalIdxs.size();++i)
            if(globalBlockDofsFilter[i])
              blocksIDs_[globalIdxs[i]]=boundaryID;
        }
      }
    }
  }
};

// Dirichlet implementation
template<typename DiscreteSpaceImp,typename CoupledMeshManagerImp>
class DirichletCondition:
  public BoundaryCondition<DiscreteSpaceImp,CoupledMeshManagerImp,int,DirichletCondition<DiscreteSpaceImp,CoupledMeshManagerImp>>
{
  public:
  typedef DiscreteSpaceImp DiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef DirichletCondition<DiscreteSpaceType,CoupledMeshManagerType> ThisType;
  typedef BoundaryCondition<DiscreteSpaceType,CoupledMeshManagerType,int,ThisType> BaseType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalAnalyticalFunctionType LocalAnalyticalFunctionType;

  friend BaseType;

  DirichletCondition(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,dirichlet)
  {}

  using BaseType::space;
  using BaseType::localBoundaryDOFs;

  template<typename... Args>
  void addBC(Args&&... args)
  {
    g_.emplace(std::forward<Args>(args)...);
  }

  void addAllBoundaryIDs(const LocalAnalyticalFunctionType& analyticalFunction)
  {
    for(const auto& boundaryID:meshmanager_.listUniqueIDs())
      addBC(boundaryID,analyticalFunction);
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"Dirichlet condition on IDs:\n";
    for(auto& mapEntry:g_)
      s<<mapEntry.first<<"\n";
  }

  private:
  using BaseType::blocksIDs_;
  using BaseType::g_;
  using BaseType::meshmanager_;

  // clear as unit rows the ones which have a Dirichlet condition (first operator is the one on the diagonal)
  // set in the RHS vector the Dirichlet component to the correct value
  template<typename TimeProviderType,typename RHSType,typename... OperatorsType>
  void setDOFs(const EntityType& entity,const TimeProviderType& timeProvider,RHSType& rhs,OperatorsType&... operators) const
  {
    constexpr std::size_t operatorsSize(sizeof...(OperatorsType));
    auto contained(isEntityContained(entity,operators.domainSpace().gridPart()...));
    auto localMatrices(std::make_tuple(operators.systemMatrix().localMatrix()...));
    auto localRHSs(getUninitializedLocalFunctions(rhs));
    Hybrid::forEach(std::make_index_sequence<operatorsSize>{},
      [&](auto i)
      {
        if(std::get<i>(contained))
        {
          std::get<i>(localMatrices).init(entity,entity);
          std::get<i>(localRHSs).init(entity);
        }
      });;
    constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
    const auto numLocalBlocks(space().basisFunctionSet(entity).size()/DiscreteSpaceType::dimRange);
    std::vector<std::size_t> globalIdxs(numLocalBlocks);
    space().blockMapper().map(entity,globalIdxs);

    std::size_t row(0);
    for(auto localIdx=decltype(numLocalBlocks){0};localIdx!=numLocalBlocks;++localIdx)
    {
      const auto& boundaryID(blocksIDs_[globalIdxs[localIdx]]);
      if(boundaryID>-1)
      {
        const auto& localBCDOFs(localBoundaryDOFs(timeProvider.time(),entity,boundaryID));
        for(auto l=decltype(localBlockSize){0};l!=localBlockSize;++l,++row)
        {
          // impose bc on first operator
          auto& firstLocalMatrix(std::get<0>(localMatrices));
          auto& firstLocalRHS(std::get<0>(localRHSs));
          firstLocalMatrix.clearRow(row);
          #if USE_SYMMETRIC_DIRICHLET
          for(auto j=decltype(firstLocalMatrix.rows()){0};j!=firstLocalMatrix.rows();++j)
          {
            // keep matrix symmetric
            firstLocalRHS[j]-=firstLocalMatrix.get(j,row)*localBCDOFs[row];
            firstLocalMatrix.set(j,row,0.0);
          }
          #endif
          firstLocalMatrix.set(row,row,1.0);
          // impose bc on others operators
          #if USE_SYMMETRIC_DIRICHLET
          Hybrid::forEach(std::make_index_sequence<operatorsSize-1>{},
            [&](auto i)
              {
                if(std::get<i+1>(contained))
                {
                  auto& currentLocalMatrix(std::get<i+1>(localMatrices));
                  auto& currentLocalRHS(std::get<i+1>(localRHSs));
                  for(auto j=decltype(currentLocalMatrix.columns()){0};j!=currentLocalMatrix.columns();++j)
                  {
                    currentLocalRHS[j]-=currentLocalMatrix.get(row,j)*localBCDOFs[row];
                    currentLocalMatrix.set(row,j,0.0);
                  }
                }
              });
          #else
          Hybrid::forEach(std::make_index_sequence<operatorsSize-1>{},
            [&](auto i)
            {
              if(std::get<i+1>(contained))
                std::get<i+1>(localMatrices).clearRow(row);
            });
          #endif
          // impose bc on RHS term
          firstLocalRHS[row]=localBCDOFs[row];
        }
      }
      else
        row+=localBlockSize;
    }
  }
};

// Free-slip implementation
template<typename DiscreteSpaceImp,typename CoupledMeshManagerImp>
class FreeSlipCondition:
  public BoundaryCondition<DiscreteSpaceImp,CoupledMeshManagerImp,ListBlockID,FreeSlipCondition<DiscreteSpaceImp,CoupledMeshManagerImp>>
{
  public:
  typedef DiscreteSpaceImp DiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef FreeSlipCondition<DiscreteSpaceType,CoupledMeshManagerType> ThisType;
  typedef BoundaryCondition<DiscreteSpaceType,CoupledMeshManagerType,ListBlockID,ThisType> BaseType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::LocalAnalyticalFunctionType LocalAnalyticalFunctionType;

  friend BaseType;

  FreeSlipCondition(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,freeslip)
  {}

  using BaseType::evaluateBoundaryFunction;
  using BaseType::space;
  using BaseType::localBoundaryDOFs;

  void addBC(int boundaryID)
  {
    for(const auto& entity:elements(meshmanager_.bulkGridPart()))
      for(const auto& intersection:intersections(meshmanager_.bulkGridPart(),entity))
        if(intersection.boundary())
          if(boundaryID==meshmanager_.intersectionID(intersection))
          {
            RangeType value(1);
            const auto normalVector(intersection.centerUnitOuterNormal());
            for(auto i=decltype(value.size()){0};i!=value.size();++i)
            {
              RangeType e(0);
              e[i]=1;
              value[i]-=std::round(std::abs(e.dot(normalVector)));
              if((value[i]!=0.0)&&(value[i]!=1.0))
                DUNE_THROW(InvalidStateException,"ERROR: free-slip function can have only 0 or 1 value!");
            }
            g_.emplace(boundaryID,LocalAnalyticalFunctionType([val=value](const DomainType& ,double ,const EntityType& ){return val;}));
            return;
          }
  }

  void addAllBoundaryIDs()
  {
    for(const auto& boundaryID:meshmanager_.listUniqueIDs())
      addBC(boundaryID);
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    const auto entity(*(meshmanager_.bulkGridPart().template begin<0>()));
    const auto x(entity.geometry().center());
    s<<"Free-slip condition on IDs:\n";
    for(auto& mapEntry:g_)
    {
      mapEntry.second.init(entity);
      mapEntry.second.initialize(0,0);
      RangeType ret;
      mapEntry.second.evaluate(x,ret);
      s<<mapEntry.first<<" ( multiplication mask "<<ret<<" )\n";
    }
  }

  private:
  using BaseType::blocksIDs_;
  using BaseType::g_;
  using BaseType::meshmanager_;

  // clear as unit rows the one which have a free-slip condition (first operator is the one on the diagonal)
  // set in the RHS vector the free-slip condition
  template<typename Foo,typename RHSType,typename... OperatorsType>
  void setDOFs(const EntityType& entity,const Foo& ,RHSType& rhs,OperatorsType&... operators) const
  {
    constexpr std::size_t operatorsSize(sizeof...(OperatorsType));
    auto contained(isEntityContained(entity,operators.domainSpace().gridPart()...));
    auto localMatrices(std::make_tuple(operators.systemMatrix().localMatrix()...));
    auto localRHSs(getUninitializedLocalFunctions(rhs));
    std::get<0>(localRHSs).init(entity);
    Hybrid::forEach(std::make_index_sequence<operatorsSize>{},
      [&](auto i)
      {
        if(std::get<i>(contained))
          std::get<i>(localMatrices).init(entity,entity);
      });
    constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
    const auto numLocalBlocks(space().basisFunctionSet(entity).size()/DiscreteSpaceType::dimRange);
    std::vector<std::size_t> globalIdxs(numLocalBlocks);
    space().blockMapper().map(entity,globalIdxs);

    std::size_t row(0);
    for(auto localIdx=decltype(numLocalBlocks){0};localIdx!=numLocalBlocks;++localIdx)
    {
      const auto& boundaryIDs(blocksIDs_[globalIdxs[localIdx]].get());
      if(boundaryIDs.front()>-1)
      {
        for(const auto& boundaryID:boundaryIDs)
        {
          const auto& localBCDOFs(localBoundaryDOFs(0.0,entity,boundaryID));
          const auto f(evaluateBoundaryFunction(typename BaseType::DomainType(0.0),0.0,entity,boundaryID));
          for(auto l=decltype(localBlockSize){0};l!=localBlockSize;++l,++row)
          {
            // impose bc on operators
            if(f[l]==0.0)
            {
              Hybrid::forEach(std::make_index_sequence<operatorsSize>{},
                [&](auto i)
                {
                  if(std::get<i>(contained))
                    std::get<i>(localMatrices).clearRow(row);
                });
              std::get<0>(localMatrices).set(row,row,1.0);
            }
            // impose bc on RHS term
            std::get<0>(localRHSs)[row]*=localBCDOFs[row];
          }
          row-=localBlockSize;
        }
      }
      row+=localBlockSize;
    }
  }
};

}
}

#endif // DUNE_FEM_BOUNDARYCONDITION_HH
