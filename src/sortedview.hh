#ifndef DUNE_FEM_SORTED_VIEW_HH
#define DUNE_FEM_SORTED_VIEW_HH

#include <array>
#include <cmath>
#include <list>
#include <vector>

#include <dune/common/fvector.hh>

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

  SortedView(const GridType& grid,double epsilon,const FieldVector<double,dimensionworld>& x0,const FieldVector<double,dimensionworld>& x1):
    grid_(grid),epsilon_(epsilon)
  {
    // compute sizes of reticulus
    const auto dx(x1-x0);
    std::array<unsigned int,3> size{1,1,1};
    for(std::size_t i=0;i!=dimensionworld;++i)
      size[i]+=static_cast<unsigned int>(ceil(dx[i]/epsilon_));
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
        idx[i]=static_cast<unsigned int>(std::round((centre[i]-x0[i])/epsilon_));
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
  const double epsilon_;
  DataType entities_;
};

}
}

#endif // DUNE_FEM_SORTED_VIEW_HH
