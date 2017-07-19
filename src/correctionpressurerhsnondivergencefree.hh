#ifndef DUNE_FEM_CORRECTIONPRESSURERHSNONDIVERGENCEFREE_HH
#define DUNE_FEM_CORRECTIONPRESSURERHSNONDIVERGENCEFREE_HH

#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/integrator.hh>

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

  // compute (fdiv,1)/(1,1)
  Integrator<CachingQuadrature<typename DiscreteFunctionType::GridPartType,0>> integrator(2*space.order()+1);
  RangeType temp(0);
  typename DiscreteFunctionType::RangeFieldType vol(0);
  for(const auto& entity:space)
  {
    vol+=std::abs(entity.geometry().volume());
    problem.pressureRHS().init(entity);
    integrator.integrateAdd(entity,problem.pressureRHS(),temp);
  }
  temp/=vol;

  // compute (fdiv-temp,\phi) and apply correction to the RHS
  LocalContribution<DiscreteFunctionType,Assembly::Add> localRHS(rhs);
  for(const auto& entity:space)
  {
    problem.pressureRHS().init(entity);
    localRHS.bind(entity);
    const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space.order()+1);
    for(const auto& qp:quadrature)
    {
      RangeType fdivValue;
      problem.pressureRHS().evaluate(qp.position(),fdivValue);
      space.basisFunctionSet(entity).evaluateAll(qp,phi);
      const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
      for(auto row=decltype(localRHS.size()){0};row!=localRHS.size();++row)
        localRHS[row]+=weight*(fdivValue-temp)*phi[row];
    }
    localRHS.unbind();
  }
}

}
}

#endif // DUNE_FEM_CORRECTIONPRESSURERHSNONDIVERGENCEFREE_HH
