#ifndef DUNE_FEM_SORTED_VIEW_HH
#define DUNE_FEM_SORTED_VIEW_HH

#include <array>
#include <cmath>
#include <list>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename SortedViewImp>
class EntityIteratorSortedView
{
  public:
  typedef SortedViewImp SortedViewType;
  typedef EntityIteratorSortedView<SortedViewType> ThisType;

  typedef typename SortedViewType::DataType DataType;

  EntityIteratorSortedView()=delete;
  EntityIteratorSortedView(const EntityIteratorSortedView& )=default;
  EntityIteratorSortedView& operator=(const EntityIteratorSortedView& )=default;

  EntityIteratorSortedView(const SortedViewType& view,bool isEnd=false):
    view_(view),datait_(isEnd?view_.data().end():view_.data().begin()),dataitend_(view_.data().end())
  {
    // set seed iterator to the first valid seed
    for(;datait_!=dataitend_;++datait_)
    {
      seedit_=datait_->begin();
      if(seedit_!=datait_->end())
        break;
    }
  }

  ThisType& operator++()
  {
    ++seedit_;
    if(seedit_==datait_->end())
      for(++datait_;datait_!=dataitend_;++datait_)
      {
        seedit_=datait_->begin();
        if(seedit_!=datait_->end())
          break;
      }
    return *this;
  }

  bool operator==(const ThisType& other) const
  {
    bool equal(false);
    if((datait_==dataitend_)&&(other.datait_==dataitend_))
      equal=true;
    else
    {
      if((datait_!=dataitend_)&&(other.datait_!=dataitend_))
      {
        if(seedit_==other.seedit_)
          equal=true;
      }
    }
    return equal;
  }

  bool operator!=(const ThisType& other) const
  {
    return !((*this)==other);
  }

  typename SortedViewType::EntityType operator*() const
  {
    return view_.grid().entity(*seedit_);
  }

  private:
  const SortedViewType& view_;
  typename DataType::const_iterator datait_;
  typename DataType::const_iterator dataitend_;
  typename DataType::value_type::const_iterator seedit_;
};

template<typename GridImp>
class SortedView
{
  public:
  typedef GridImp GridType;
  typedef SortedView<GridType> ThisType;

  static constexpr auto dimensionworld=GridType::dimensionworld;
  static constexpr auto dimension=GridType::dimension;

  typedef typename GridType::template Codim<0>::Entity EntityType;
  typedef typename EntityType::EntitySeed EntitySeedType;
  typedef std::vector<std::list<EntitySeedType>> DataType;
  typedef EntityIteratorSortedView<ThisType> Iterator;

  template<typename BulkBoundingBoxType>
  SortedView(const GridType& grid,const BulkBoundingBoxType& bulkBoundingBox):
    grid_(grid)
  {
    // epsilon is the length of the edge/face of the equilateral triangle/tetrahedron which has volume equal to the total volume of
    // the bounding box divivided the actual number of entities therefore
    // epsilon = (bounding_box_volume / num_entities)^(1/dimensionworld)*k where k is 1.5 if dimensionworld==2 and 2 if dimensionworld==3
    const auto dx(bulkBoundingBox.second-bulkBoundingBox.first);
    double boundingBoxVolume(1.0);
    for(std::size_t i=0;i!=dimensionworld;++i)
      boundingBoxVolume*=dx[i];
    double epsilon(std::pow(boundingBoxVolume/static_cast<double>(grid.size(0)),1.0/static_cast<double>(dimensionworld)));
    epsilon*=(dimensionworld==2?1.5:2.0);
    // compute sizes of reticulus
    std::array<unsigned int,3> size{1,1,1};
    for(std::size_t i=0;i!=dimensionworld;++i)
      size[i]+=static_cast<unsigned int>(std::ceil(dx[i]/epsilon));
    // resize vector of entities
    entities_.resize(size[0]*size[1]*size[2]);
    // fill entities
    const auto& leafGridView(grid.leafGridView());
    for(const auto& entity:elements(leafGridView))
    {
      // compute barycenter
      const auto centre(entity.geometry().center());
      // compute position inside entitites
      std::array<unsigned int,3> idx{0,0,0};
      for(std::size_t i=0;i!=dimensionworld;++i)
        idx[i]=static_cast<unsigned int>(std::round((centre[i]-bulkBoundingBox.first[i])/epsilon));
      auto position(idx[2]*size[0]*size[1]);
      if(idx[2]%2==0)
      {
        position+=((idx[1]%2==0)?(idx[1]*size[0]+idx[0]):((idx[1]+1)*size[0]-idx[0]-1));
      }
      else
      {
        position+=size[0]*size[1]-1;
        position-=((idx[1]%2==0)?(idx[1]*size[0]+idx[0]):((idx[1]+1)*size[0]-idx[0]-1));
      }
      // add entity seed into entities
      entities_[position].emplace_back(std::move(entity.seed()));
    }
  }

  const GridType& grid() const
  {
    return grid_;
  }

  const DataType& data() const
  {
    return entities_;
  }

  template<int codim>
  Iterator begin() const
  {
    static_assert(codim==0,"EntityIteratorSortedView exists only for codim=0!");
    return Iterator(*this);
  }

  template<int codim>
  Iterator end() const
  {
    static_assert(codim==0,"EntityIteratorSortedView exists only for codim=0!");
    return Iterator(*this,true);
  }

  private:
  const GridType& grid_;
  DataType entities_;
};

}
}

#endif // DUNE_FEM_SORTED_VIEW_HH
