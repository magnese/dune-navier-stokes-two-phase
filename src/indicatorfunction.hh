#ifndef DUNE_INDICATORFUNCTION_HH
#define DUNE_INDICATORFUNCTION_HH

#include <vector>

namespace Dune
{

// indicator function of the inner part of the bulk grid
template<class BulkGridImp>
class IndicatorFunction
{
  public:
  typedef BulkGridImp BulkGridType;
  typedef IndicatorFunction<BulkGridType> ThisType;

  // constructor
  explicit IndicatorFunction(const BulkGridType& bulkGrid,const std::vector<int>& elementIDs):
    bulkgrid_(bulkGrid),elementids_(elementIDs)
  {}

  IndicatorFunction(const ThisType& )=delete;

  // check if an entity is inner or outer
  template<class EntityType>
  inline bool isInner(const EntityType& entity) const
  {
    const auto idx(bulkgrid_.levelIndexSet(0).index(entity));
    return elementids_[idx]==1?true:false;
  }

  // return 1.0 if inner or 0.0 if outer
  template<class EntityType>
  inline double operator()(const EntityType& entity) const
  {
    const auto idx(bulkgrid_.levelIndexSet(0).index(entity));
    return elementids_[idx]==1?1.0:0.0;
  }

  private:
  const BulkGridType& bulkgrid_;
  const std::vector<int>& elementids_;
};

}

#endif // DUNE_INDICATORFUNCTION_HH
