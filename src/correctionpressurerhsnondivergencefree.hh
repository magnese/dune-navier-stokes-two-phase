#ifndef DUNE_FEM_CORRECTIONPRESSURERHSNONDIVERGENCEFREE_HH
#define DUNE_FEM_CORRECTIONPRESSURERHSNONDIVERGENCEFREE_HH

#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <cmath>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionType,typename ProblemType,typename TimeProviderType>
void correctionPressureRHSNonDivergenceFree(DiscreteFunctionType& rhs,ProblemType& problem,const TimeProviderType& timeProvider)
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  const auto& space(rhs.space());
  typedef typename DiscreteFunctionType::RangeType RangeType;
  std::vector<RangeType> phi(space.maxNumDofs());
  problem.pressureRHS().initialize(timeProvider.time(),timeProvider.time());

  // compute domain volume
  double vol(0.0);
  for(const auto& entity:space)
    vol+=std::abs(entity.geometry().volume());

  // perform a grid walkthrough and apply correction to the RHS
  LocalContribution<DiscreteFunctionType,Assembly::Add> localRHS(rhs);
  for(const auto& entity:space)
  {
    problem.pressureRHS().init(entity);
    localRHS.bind(entity);
    const auto& baseSet(space.basisFunctionSet(entity));
    RangeType coeff(0);

    const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space.order()+1);
    for(const auto& qp:quadrature)
    {
      baseSet.evaluateAll(qp,phi);
      const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
      for(auto row=decltype(localRHS.size()){0};row!=localRHS.size();++row)
        coeff+=phi[row]*weight;
    }
    coeff/=vol;
    for(const auto& qp:quadrature)
    {
      RangeType fdivValue;
      problem.pressureRHS().evaluate(qp.position(),fdivValue);
      baseSet.evaluateAll(qp,phi);
      const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
      for(auto row=decltype(localRHS.size()){0};row!=localRHS.size();++row)
      {
        auto value(fdivValue*(phi[row]-coeff));
        value*=weight;
        localRHS[row]+=value;
      }
    }
    localRHS.unbind();
  }
}

}
}

#endif // DUNE_FEM_CORRECTIONPRESSURERHSNONDIVERGENCEFREE_HH
