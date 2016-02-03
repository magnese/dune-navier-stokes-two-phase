#ifndef DUNE_INDICATORFUNCTION_HH
#define DUNE_INDICATORFUNCTION_HH

#include <vector>

namespace Dune
{

// indicator function of the inner part of the bulk grid
template<typename BulkGridImp>
class IndicatorFunction
{
  public:
  typedef BulkGridImp BulkGridType;
  typedef IndicatorFunction<BulkGridType> ThisType;

  typedef typename BulkGridType::template Codim<0>::Entity BulkEntityType;
  typedef typename BulkGridType::HostGrid::template Codim<0>::Entity HostEntityType;

  // constructor
  explicit IndicatorFunction(const BulkGridType& bulkGrid,const std::vector<int>& elementIDs):
    bulkgrid_(bulkGrid),elementids_(elementIDs)
  {}

  IndicatorFunction(const ThisType& )=delete;

  // check if the host entity is inner or outer
  bool isInner(const HostEntityType& entity) const
  {
    const auto idx(bulkgrid_.hostGrid().levelIndexSet(0).index(entity));
    return elementids_[idx]==1?true:false;
  }

  // check if the bulk entity is inner or outer
  bool isInner(const BulkEntityType& entity) const
  {
    const auto idx(bulkgrid_.levelIndexSet(0).index(entity));
    return elementids_[idx]==1?true:false;
  }

  // return 1.0 if the entity is inner or 0.0 if outer
  template<typename EntityType>
  double operator()(const EntityType& entity) const
  {
    return static_cast<double>(isInner(entity));
  }

  private:
  const BulkGridType& bulkgrid_;
  const std::vector<int>& elementids_;
};

}

#endif // DUNE_INDICATORFUNCTION_HH
