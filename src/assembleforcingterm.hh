#ifndef DUNE_FEM_ASSEMBLEFORCINGTERM_HH
#define DUNE_FEM_ASSEMBLEFORCINGTERM_HH

#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionType,typename FluidStateType,typename ProblemType,typename TimeProviderType>
void assembleForcingTerm(DiscreteFunctionType& rhs,const FluidStateType& fluidState,ProblemType& problem,
                         const TimeProviderType& timeProvider)
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  const auto& space(rhs.space());
  std::vector<typename DiscreteFunctionType::RangeType> phi(space.maxNumDofs());
  problem.velocityRHS().initialize(timeProvider.time(),timeProvider.time());

  // perform a grid walkthrough and assemble the forcing term
  LocalContribution<DiscreteFunctionType,Assembly::Add> localRHS(rhs);
  for(const auto& entity:space)
  {
    problem.velocityRHS().init(entity);
    localRHS.bind(entity);
    const auto& baseSet(space.basisFunctionSet(entity));

    const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space.order()+1);
    for(const auto& qp:quadrature)
    {
      typename FluidStateType::VelocityDiscreteFunctionType::RangeType fValue;
      problem.velocityRHS().evaluate(qp.position(),fValue);
      baseSet.evaluateAll(qp,phi);
      const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

      const auto localSize(localRHS.size());
      for(auto row=decltype(localSize){0};row!=localSize;++row)
        localRHS[row]+=weight*(fValue*phi[row]);
    }
    localRHS.unbind();
  }
}

}
}

#endif // DUNE_FEM_ASSEMBLEFORCINGTERM_HH
