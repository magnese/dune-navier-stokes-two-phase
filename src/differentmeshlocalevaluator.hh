#ifndef DUNE_FEM_DIFFERENTMESHLOCALEVALUATOR_HH
#define DUNE_FEM_DIFFERENTMESHLOCALEVALUATOR_HH

#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/fem/common/coordinate.hh>
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
    df_(df),oldentity_(*(df_.space().begin())),searchiterations_(0)
  {}

  template<class PointType>
  void evaluate(const PointType& x,RangeType& ret) const
  {
    const auto xGlobal(entity().geometry().global(coordinate(x)));
    if(useBarycentricEntitySearch)
      searchiterations_+=barycentricEntitySearch(df_.gridPart(),std::move(oldentity_),xGlobal);
    else
    {
      EntitySearch<GridPartType> entitySearch(df_.gridPart());
      oldentity_=entitySearch(xGlobal);
    }
    const auto& localDF(df_.localFunction(oldentity_));
    localDF.evaluate(oldentity_.geometry().local(xGlobal),ret);
  }

  void init(const EntityType& entity)
  {
    entity_=&entity;
  }

  const EntityType& entity() const
  {
    return *entity_;
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
};

}
}

#endif // DUNE_FEM_DIFFERENTMESHLOCALEVALUATOR_HH
