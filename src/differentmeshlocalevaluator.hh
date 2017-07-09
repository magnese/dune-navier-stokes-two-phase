#ifndef DUNE_FEM_DIFFERENTMESHLOCALEVALUATOR_HH
#define DUNE_FEM_DIFFERENTMESHLOCALEVALUATOR_HH

#include <dune/common/exceptions.hh>
#include <dune/fem/common/coordinate.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>

#include "barycentricentitysearch.hh"

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionImpl,bool useBarycentricEntitySearch=true>
class DifferentMeshLocalEvaluator
{
  public:
  typedef DiscreteFunctionImpl DiscreteFunctionType;
  typedef DifferentMeshLocalEvaluator<DiscreteFunctionType,useBarycentricEntitySearch> ThisType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  DifferentMeshLocalEvaluator(const DiscreteFunctionType& df):
    df_(df),oldentity_(*(df_.space().begin())),searchiterations_(0),localdf_(df_)
  {}

  template<class PointType>
  void evaluate(const PointType& x,RangeType& ret) const
  {
    const auto xGlobal(entity().geometry().global(coordinate(x)));
    if(useBarycentricEntitySearch)
      searchiterations_+=barycentricEntitySearch(df_.gridPart(),oldentity_,xGlobal);
    else
    {
      EntitySearch<GridPartType> entitySearch(df_.gridPart());
      oldentity_=entitySearch(xGlobal);
    }
    localdf_.init(oldentity_);
    localdf_.evaluate(oldentity_.geometry().local(xGlobal),ret);
  }

  unsigned int order() const
  {
    return df_.space().order();
  }

  void init(const EntityType& entity)
  {
    entity_=&entity;
  }

  void init(const EntityType& entiy,const EntityType& oldEntity)
  {
    init(entiy);
    oldentity_=oldEntity;
  }

  const EntityType& entity() const
  {
    return *entity_;
  }

  const EntityType& oldEntity() const
  {
    return oldentity_;
  }

  unsigned int searchIterations() const
  {
    return searchiterations_;
  }

  private:
  const DiscreteFunctionType& df_;
  mutable EntityType oldentity_;
  mutable unsigned int searchiterations_;
  const EntityType* entity_;
  mutable ConstLocalDiscreteFunction<DiscreteFunctionType> localdf_;
};

}
}

#endif // DUNE_FEM_DIFFERENTMESHLOCALEVALUATOR_HH
