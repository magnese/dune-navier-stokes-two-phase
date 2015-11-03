#ifndef DUNE_BOUNDARYCONDITION_HH
#define DUNE_BOUNDARYCONDITION_HH

#include <vector>
#include <map>
#include <utility>
#include <functional>
#include <tuple>
#include <list>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>

namespace Dune
{
namespace Fem
{

// struct to store block IDs when overlapping IDs can NOT be handled using only one of the IDs involved
class ListBlockID
{
  public:
  typedef std::list<int> value_type;

  ListBlockID(const int& val=-1):
    values_(1,val)
  {}
  inline void clear()
  {
    values_.clear();
  }
  inline void operator=(const int& value)
  {
    if(values_.front()==-1)
      values_.front()=value;
    else
      values_.push_back(value);
  }
  inline const value_type& get() const
  {
    return values_;
  }
  private:
  value_type values_;
};

// BC  type
enum BCEnumType{dirichlet,neumann,robin,freeslip};

// BC interface
template<typename DSImp,typename RSImp,typename CMMImp,typename BIDImp,typename BCImp>
class BoundaryCondition
{
  public:
  typedef DSImp DomainSpaceType;
  typedef RSImp RangeSpaceType;
  typedef CMMImp CoupledMeshManagerType;
  typedef BCImp BCImplementation;
  typedef BIDImp BlockIDType;
  typedef BoundaryCondition<DomainSpaceType,RangeSpaceType,CoupledMeshManagerType,BlockIDType,BCImplementation> ThisType;

  typedef typename DomainSpaceType::GridPartType GridPartType;
  typedef typename DomainSpaceType::DomainType DomainType;
  typedef typename DomainSpaceType::EntityType EntityType;
  typedef typename GridPartType::IntersectionIteratorType::Intersection IntersectionType;
  typedef typename RangeSpaceType::RangeType RangeType;

  typedef std::function<RangeType(const DomainType&,const double&,const EntityType&)> FunctionType;
  typedef LocalAnalyticalFunctionBinder<DomainSpaceType> LocalAnalyticalFunctionType;
  typedef std::map<int,LocalAnalyticalFunctionType> FunctionMapType;
  typedef LocalFunctionAdapter<LocalAnalyticalFunctionType> AdaptedDiscreteFunctionType;
  typedef std::map<int,AdaptedDiscreteFunctionType> AdaptedFunctionMapType;

  inline void addBC(int&& boundaryID,FunctionType&& g)
  {
    g_.emplace(boundaryID,g);
  }

  inline void addBC(int&& boundaryID,const FunctionType& g)
  {
    g_.emplace(boundaryID,g);
  }

  inline const DomainSpaceType& domainSpace() const
  {
    return *domainspace_;
  }

  inline const RangeSpaceType& rangeSpace() const
  {
    return *rangespace_;
  }

  template<typename... Args>
  inline void applyToOperator(Args&&... args)
  {
    updateDOFs();
    for(const auto entity:domainSpace())
      asImp().setDOFsMatrix(entity,args...);
  }

  template<typename DiscreteFunctionType,typename... Args>
  inline void applyToRHS(DiscreteFunctionType& w,Args&&... args)
  {
    updateDOFs();
    DiscreteFunctionType g("g",w.space());
    for(const auto entity:domainSpace())
      asImp().setDOFsRHS(entity,w,g,args...);
  }

  inline RangeType evaluateBoundaryFunction(const DomainType& x,const double& t,const EntityType& entity,const int& boundaryID) const
  {
    auto& g(g_.find(boundaryID)->second);
    g.init(entity,t);
    RangeType ret;
    g.evaluate(x,ret);
    return ret;
  }

  inline RangeType evaluateBoundaryFunction(const DomainType& x,const double& t,const IntersectionType& intersection) const
  {
    const auto& boundaryID(meshmanager_.boundaryIDs()[intersection.boundarySegmentIndex()]);
    return evaluateBoundaryFunction(x,t,intersection.inside(),boundaryID);
  }

  template<typename DiscreteFunctionType>
  void localInterpolateBoundaryFunction(const double& t,const EntityType& entity,const int& boundaryID,DiscreteFunctionType& df) const
  {
    auto gIt(gadapted_.find(boundaryID));
    if(gIt==gadapted_.end())
      DUNE_THROW(RangeError,"boundary ID not found in BC");
    auto& gAdapted(gIt->second);
    gAdapted.localFunctionImpl().initialize(t);
    const auto interpolation(df.space().interpolation(entity));
    typedef typename AdaptedDiscreteFunctionType::RangeFieldType RangeFieldType;
    std::vector<RangeFieldType> localDOFs(df.space().basisFunctionSet(entity).size());
    interpolation(gAdapted.localFunction(entity),localDOFs);
    df.setLocalDofs(entity,localDOFs);
  }

  template<typename DiscreteFunctionType>
  inline void localInterpolateBoundaryFunction(const double& t,const IntersectionType& intersection,DiscreteFunctionType& df) const
  {
    const auto& boundaryID(meshmanager_.boundaryIDs()[intersection.boundarySegmentIndex()]);
    localInterpolateBoundaryFunction(t,intersection.inside(),boundaryID,df);
  }

  BoundaryCondition()=delete;
  BoundaryCondition(const ThisType&)=delete;
  ThisType& operator=(const ThisType&)=delete;
  BoundaryCondition(ThisType&&)=default;
  ThisType& operator=(ThisType&&)=default;

  protected:
  CoupledMeshManagerType& meshmanager_;
  std::unique_ptr<GridPartType> gridpart_;
  std::unique_ptr<DomainSpaceType> domainspace_;
  std::unique_ptr<RangeSpaceType> rangespace_;
  const BCEnumType bctype_;
  mutable FunctionMapType g_;
  mutable AdaptedFunctionMapType gadapted_;
  std::vector<BlockIDType> blocksIDs_;
  unsigned int sequence_;

  BoundaryCondition(CoupledMeshManagerType& meshManager,BCEnumType bctype):
    meshmanager_(meshManager),bctype_(bctype),g_(),gadapted_(),blocksIDs_(0),sequence_(0)
  {}

  inline BCImplementation& asImp()
  {
    return static_cast<BCImplementation&>(*this);
  }

  inline const BCImplementation& asImp() const
  {
    return static_cast<const BCImplementation&>(*this);
  }

  // update the blocks IDs when the mesh changes
  void updateDOFs()
  {
    // check if the mesh is changed
    if(sequence_!=meshmanager_.sequence())
    {
      // reset space pointers (to avoid segmentation fault)
      rangespace_.reset();
      domainspace_.reset();
      gridpart_.reset();
      // create grid part and spaces
      gridpart_=std::unique_ptr<GridPartType>(new GridPartType(meshmanager_.bulkGrid()));
      domainspace_=std::unique_ptr<DomainSpaceType>(new DomainSpaceType(*gridpart_));
      rangespace_=std::unique_ptr<RangeSpaceType>(new RangeSpaceType(*gridpart_));
      // set blocksIDs
      blocksIDs_.clear();
      blocksIDs_.resize(domainSpace().blockMapper().size(),-1);
      for(const auto entity:domainSpace())
      {
        if(entity.hasBoundaryIntersections())
          setBlocksIDs(entity);
      }
      // create adapted discrete boundary function
      gadapted_.clear();
      for(auto& g:g_)
        gadapted_.emplace(g.first,AdaptedDiscreteFunctionType("adapted function",g.second,*gridpart_));
      // update sequence number
      sequence_=meshmanager_.sequence();
    }
  }

  // set the correct IDs for each block
  void setBlocksIDs(const EntityType& entity)
  {
    const auto& blockMapper(domainSpace().blockMapper());
    std::vector<std::size_t> globalIdxs(blockMapper.numDofs(entity));
    blockMapper.map(entity,globalIdxs);
    std::vector<bool> globalBlockDofsFilter(blockMapper.numDofs(entity));

    for(const auto intersection:intersections(static_cast<typename GridPartType::GridViewType>(*gridpart_),entity))
    {
      if(intersection.boundary())
      {
        const auto& boundaryID(meshmanager_.boundaryIDs()[intersection.boundarySegmentIndex()]);
        if(g_.find(boundaryID)!=g_.end())
        {
          const auto faceLocalIdx(intersection.indexInInside());
          blockMapper.onSubEntity(entity,faceLocalIdx,1,globalBlockDofsFilter);
          for(auto i=0;i!=globalIdxs.size();++i)
            if(globalBlockDofsFilter[i])
              blocksIDs_[globalIdxs[i]]=boundaryID;
        }
      }
    }
  }

  class ClearRows
  {
    public:
    template<typename T>
    ClearRows(T& tup,const std::size_t& row)
    {
      constexpr auto size(std::tuple_size<T>::value);
      Caller<T,size> caller(tup,row);
    }

    private:
    template<typename T,std::size_t n>
    struct Caller
    {
      inline Caller(T& tup,const std::size_t& row)
      {
        std::get<n-1>(tup).clearRow(row);
        Caller<T,n-1>(tup,row);
      }
    };
    template<typename T>
    struct Caller<T,0>
    {
      inline Caller(T& ,const std::size_t& )
      {}
    };
  };
};

// Dirichlet implementation
template<typename DSImp,typename RSImp,typename CMMImp>
class DirichletCondition:public BoundaryCondition<DSImp,RSImp,CMMImp,int,DirichletCondition<DSImp,RSImp,CMMImp>>
{
  public:
  typedef DSImp DomainSpaceType;
  typedef RSImp RangeSpaceType;
  typedef CMMImp CoupledMeshManagerType;
  typedef DirichletCondition<DomainSpaceType,RangeSpaceType,CoupledMeshManagerType> ThisType;
  typedef BoundaryCondition<DomainSpaceType,RangeSpaceType,CoupledMeshManagerType,int,ThisType> BaseType;
  typedef typename BaseType::EntityType EntityType;

  friend class BoundaryCondition<DomainSpaceType,RangeSpaceType,CoupledMeshManagerType,int,ThisType>;

  DirichletCondition(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,dirichlet)
  {}

  using BaseType::evaluateBoundaryFunction;
  using BaseType::domainSpace;
  using BaseType::localInterpolateBoundaryFunction;

  private:
  using BaseType::blocksIDs_;

  // clear as unit rows the ones which have a Dirichlet condition (first operator is the one on the diagonal)
  template<typename ...OperatorsType>
  void setDOFsMatrix(const EntityType& entity,OperatorsType&... operators) const
  {
    typedef std::tuple<typename OperatorsType::LinearOperatorType::LocalMatrixType...> LocalMatricesType;
    LocalMatricesType localMatrices(operators.systemMatrix().localMatrix(entity,entity)...);
    const auto localBlockSize(DomainSpaceType::localBlockSize);
    const auto numLocalBlocks(std::get<0>(localMatrices).rows()/localBlockSize);
    std::vector<std::size_t> globalIdxs(numLocalBlocks);
    domainSpace().blockMapper().map(entity,globalIdxs);

    std::size_t row(0);
    for(auto localIdx=0;localIdx!=numLocalBlocks;++localIdx)
    {
      const auto& boundaryID(blocksIDs_[globalIdxs[localIdx]]);
      if(boundaryID>-1)
      {
        for(auto l=0;l!=localBlockSize;++l,++row)
        {
          typename BaseType::ClearRows clearRows(localMatrices,row);
          std::get<0>(localMatrices).set(row,row,1.0);
        }
      }
      else
        row+=localBlockSize;
    }
  }

  // set in the RHS vector the Dirichlet component to the correct value
  template<typename DiscreteFunctionType,typename TimeProviderType>
  void setDOFsRHS(const EntityType& entity,DiscreteFunctionType& w,DiscreteFunctionType& g,
                  const TimeProviderType& timeProvider) const
  {
    auto wLocal(w.localFunction(entity));
    auto gLocal(g.localFunction(entity));
    const auto numLocalBlocks(wLocal.numScalarDofs());
    const auto localBlockSize(DiscreteFunctionType::DiscreteFunctionSpaceType::localBlockSize);
    std::vector<std::size_t> globalIdxs(numLocalBlocks);
    domainSpace().blockMapper().map(entity,globalIdxs);

    std::size_t row(0);
    for(auto localIdx=0;localIdx!=numLocalBlocks;++localIdx)
    {
      const auto& boundaryID(blocksIDs_[globalIdxs[localIdx]]);
      if(boundaryID>-1)
      {
        localInterpolateBoundaryFunction(timeProvider.time(),entity,boundaryID,g);
        for(auto l=0;l!=localBlockSize;++l,++row)
          wLocal[row]=gLocal[row];
      }
      else
        row+=localBlockSize;
    }
  }
};

// Free-slip implementation
template<typename DSImp,typename RSImp,typename CMMImp>
class FreeSlipCondition:public BoundaryCondition<DSImp,RSImp,CMMImp,ListBlockID,FreeSlipCondition<DSImp,RSImp,CMMImp>>
{
  public:
  typedef DSImp DomainSpaceType;
  typedef RSImp RangeSpaceType;
  typedef CMMImp CoupledMeshManagerType;
  typedef FreeSlipCondition<DomainSpaceType,RangeSpaceType,CoupledMeshManagerType> ThisType;
  typedef BoundaryCondition<DomainSpaceType,RangeSpaceType,CoupledMeshManagerType,ListBlockID,ThisType> BaseType;
  typedef typename BaseType::EntityType EntityType;

  friend class BoundaryCondition<DomainSpaceType,RangeSpaceType,CoupledMeshManagerType,ListBlockID,ThisType>;

  FreeSlipCondition(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,freeslip)
  {}

  using BaseType::evaluateBoundaryFunction;
  using BaseType::domainSpace;
  using BaseType::localInterpolateBoundaryFunction;

  private:
  using BaseType::blocksIDs_;

  // clear as unit rows the one which have a free-slip condition (first operator is the one on the diagonal)
  template<typename... OperatorsType>
  void setDOFsMatrix(const EntityType& entity,OperatorsType&... operators) const
  {
    typedef std::tuple<typename OperatorsType::LinearOperatorType::LocalMatrixType...> LocalMatricesType;
    LocalMatricesType localMatrices(operators.systemMatrix().localMatrix(entity,entity)...);
    const auto localBlockSize(DomainSpaceType::localBlockSize);
    const auto numLocalBlocks(std::get<0>(localMatrices).rows()/localBlockSize);
    std::vector<std::size_t> globalIdxs(numLocalBlocks);
    domainSpace().blockMapper().map(entity,globalIdxs);

    std::size_t row(0);
    for(auto localIdx=0;localIdx!=numLocalBlocks;++localIdx)
    {
      const auto& boundaryIDs(blocksIDs_[globalIdxs[localIdx]].get());
      if(boundaryIDs.front()>-1)
      {
        for(const auto& boundaryID:boundaryIDs)
        {
          const auto f(evaluateBoundaryFunction(typename BaseType::DomainType(0.0),0.0,entity,boundaryID));
          for(auto l=0;l!=localBlockSize;++l,++row)
          {
            if(f[l]==0.0)
            {
              typename BaseType::ClearRows clearRows(localMatrices,row);
              std::get<0>(localMatrices).set(row,row,1.0);
            }
          }
          row-=localBlockSize;
        }
      }
      row+=localBlockSize;
    }
  }

  // set in the RHS vector the free-slip condition
  template<typename DiscreteFunctionType,typename... Args>
  void setDOFsRHS(const EntityType& entity,DiscreteFunctionType& w,DiscreteFunctionType& g,Args&&... ) const
  {
    auto wLocal(w.localFunction(entity));
    auto gLocal(g.localFunction(entity));
    const auto numLocalBlocks(wLocal.numScalarDofs());
    const auto localBlockSize(DiscreteFunctionType::DiscreteFunctionSpaceType::localBlockSize);
    std::vector<std::size_t> globalIdxs(numLocalBlocks);
    domainSpace().blockMapper().map(entity,globalIdxs);

    std::size_t row(0);
    for(auto localIdx=0;localIdx!=numLocalBlocks;++localIdx)
    {
      const auto& boundaryIDs(blocksIDs_[globalIdxs[localIdx]].get());
      if(boundaryIDs.front()>-1)
      {
        for(const auto& boundaryID:boundaryIDs)
        {
          localInterpolateBoundaryFunction(0.0,entity,boundaryID,g);
          for(auto l=0;l!=localBlockSize;++l,++row)
            wLocal[row]*=gLocal[row];
          row-=localBlockSize;
        }
      }
      row+=localBlockSize;
    }
  }
};

}
}

#endif // DUNE_BOUNDARYCONDITION_HH
