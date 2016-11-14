#ifndef DUNE_FEM_BOUNDARYCONDITION_HH
#define DUNE_FEM_BOUNDARYCONDITION_HH

#include <functional>
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

  typedef std::function<RangeType(const DomainType&,double,const EntityType&)> FunctionType;
  typedef LocalAnalyticalFunctionBinder<DiscreteSpaceType> LocalAnalyticalFunctionType;
  typedef std::map<int,LocalAnalyticalFunctionType> FunctionMapType;
  typedef LocalFunctionAdapter<LocalAnalyticalFunctionType> AdaptedDiscreteFunctionType;
  typedef std::map<int,AdaptedDiscreteFunctionType> AdaptedFunctionMapType;

  typedef DynamicVector<typename AdaptedDiscreteFunctionType::RangeFieldType> LocalBoundaryDOFsType;

  void addBC(int boundaryID,const FunctionType& g)
  {
    g_.emplace(boundaryID,g);
  }

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
    auto& g(g_.find(boundaryID)->second);
    g.init(entity,t);
    RangeType ret;
    g.evaluate(x,ret);
    return ret;
  }

  RangeType evaluateBoundaryFunction(const DomainType& x,double t,const IntersectionType& intersection) const
  {
    const auto& boundaryID(meshmanager_.boundaryIDs()[intersection.boundarySegmentIndex()]);
    return evaluateBoundaryFunction(x,t,intersection.inside(),boundaryID);
  }

  LocalBoundaryDOFsType localBoundaryDOFs(double t,const EntityType& entity,int boundaryID) const
  {
    auto gIt(gadapted_.find(boundaryID));
    if(gIt==gadapted_.end())
      DUNE_THROW(RangeError,"boundary ID not found in BC");
    auto& gAdapted(gIt->second);
    gAdapted.localFunctionImpl().initialize(t);
    const auto interpolation(space().interpolation(entity));
    LocalBoundaryDOFsType localDOFs(space().basisFunctionSet(entity).size());
    interpolation(gAdapted.localFunction(entity),localDOFs);
    return localDOFs;
  }

  LocalBoundaryDOFsType localBoundaryDOFs(double t,const IntersectionType& intersection) const
  {
    const auto& boundaryID(meshmanager_.boundaryIDs()[intersection.boundarySegmentIndex()]);
    return localBoundaryDOFs(t,intersection.inside(),boundaryID);
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
        gadapted_.emplace(g.first,AdaptedDiscreteFunctionType("adapted function",g.second,meshmanager_.bulkGridPart()));
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
        const auto& boundaryID(meshmanager_.boundaryIDs()[intersection.boundarySegmentIndex()]);
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

  friend BaseType;

  DirichletCondition(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,dirichlet)
  {}

  using BaseType::evaluateBoundaryFunction;
  using BaseType::space;
  using BaseType::localBoundaryDOFs;

  private:
  using BaseType::blocksIDs_;

  // clear as unit rows the ones which have a Dirichlet condition (first operator is the one on the diagonal)
  // set in the RHS vector the Dirichlet component to the correct value
  template<typename TimeProviderType,typename RHSType,typename FirstOperatorType,typename... OperatorsType>
  void setDOFs(const EntityType& entity,const TimeProviderType& timeProvider,RHSType& rhs,FirstOperatorType& firstOperator,
               OperatorsType&... operators) const
  {
    typedef typename FirstOperatorType::LinearOperatorType::LocalMatrixType FirstLocalMatrixType;
    FirstLocalMatrixType firstLocalMatrix(firstOperator.systemMatrix().localMatrix(entity,entity));
    typedef std::tuple<typename OperatorsType::LinearOperatorType::LocalMatrixType...> LocalMatricesType;
    LocalMatricesType localMatrices(operators.systemMatrix().localMatrix(entity,entity)...);
    auto rhsLocal(rhs.localFunction(entity));
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
          firstLocalMatrix.clearRow(row);
          #if USE_SYMMETRIC_DIRICHLET
          for(auto j=decltype(firstLocalMatrix.rows()){0};j!=firstLocalMatrix.rows();++j)
          {
            // keep matrix symmetric
            rhsLocal[j]-=firstLocalMatrix.get(j,row)*localBCDOFs[row];
            firstLocalMatrix.set(j,row,0.0);
          }
          #endif
          firstLocalMatrix.set(row,row,1.0);
          // impose bc on others operators
          #if USE_SYMMETRIC_DIRICHLET
          std::size_t offset(firstLocalMatrix.rows());
          Hybrid::forEach(std::make_index_sequence<std::tuple_size<LocalMatricesType>::value>{},
            [&](auto i)
              {
                auto& entry(std::get<i>(localMatrices));
                for(auto j=decltype(entry.columns()){0};j!=entry.columns();++j)
                {
                  rhsLocal[j+offset]-=entry.get(row,j)*localBCDOFs[row];
                  entry.set(row,j,0.0);
                }
                offset+=entry.columns();
              });
          #else
          Hybrid::forEach(std::make_index_sequence<std::tuple_size<LocalMatricesType>::value>{},
            [&](auto i){std::get<i>(localMatrices).clearRow(row);});
          #endif
          // impose bc on RHS term
          rhsLocal[row]=localBCDOFs[row];
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

  friend BaseType;

  FreeSlipCondition(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,freeslip)
  {}

  using BaseType::evaluateBoundaryFunction;
  using BaseType::space;
  using BaseType::localBoundaryDOFs;

  private:
  using BaseType::blocksIDs_;

  // clear as unit rows the one which have a free-slip condition (first operator is the one on the diagonal)
  // set in the RHS vector the free-slip condition
  template<typename Foo,typename RHSType,typename... OperatorsType>
  void setDOFs(const EntityType& entity,const Foo& ,RHSType& rhs,OperatorsType&... operators) const
  {
    typedef std::tuple<typename OperatorsType::LinearOperatorType::LocalMatrixType...> LocalMatricesType;
    LocalMatricesType localMatrices(operators.systemMatrix().localMatrix(entity,entity)...);
    auto rhsLocal(rhs.localFunction(entity));
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
            if(f[l]==0.0)
            {
              Hybrid::forEach(std::make_index_sequence<std::tuple_size<LocalMatricesType>::value>{},
                [&](auto i){std::get<i>(localMatrices).clearRow(row);});
              std::get<0>(localMatrices).set(row,row,1.0);
            }
            rhsLocal[row]*=localBCDOFs[row];
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
