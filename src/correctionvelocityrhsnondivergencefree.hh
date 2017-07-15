#ifndef DUNE_FEM_CORRECTIONVELOCITYRHSNONDIVERGENCEFREE_HH
#define DUNE_FEM_CORRECTIONVELOCITYRHSNONDIVERGENCEFREE_HH

#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionType,typename ProblemType,typename TimeProviderType>
void correctionVelocityRHSNonDivergenceFree(DiscreteFunctionType& rhs,ProblemType& problem,const TimeProviderType& timeProvider)
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  const auto& space(rhs.space());
  typedef typename DiscreteFunctionType::RangeType RangeType;
  std::vector<RangeType> phi(space.maxNumDofs());
  problem.pressureRHS().initialize(timeProvider.time(),timeProvider.time());
  problem.velocitySolution().initialize(timeProvider.time(),timeProvider.time());

  // perform a grid walkthrough and apply correction to the RHS
  LocalContribution<DiscreteFunctionType,Assembly::Add> localRHS(rhs);
  for(const auto& entity:space)
  {
    problem.pressureRHS().init(entity);
    problem.velocitySolution().init(entity);
    localRHS.bind(entity);
    const auto& baseSet(space.basisFunctionSet(entity));

    const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space.order()+1);
    for(const auto& qp:quadrature)
    {
      typename ProblemType::PressureRangeType fdivValue;
      problem.pressureRHS().evaluate(qp.position(),fdivValue);
      RangeType velocityValue;
      problem.velocitySolution().evaluate(qp.position(),velocityValue);
      baseSet.evaluateAll(qp,phi);
      const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
      // compute rho
      const auto globalPoint(entity.geometry().global(qp.position()));
      const auto rho(globalPoint.two_norm()>problem.exactRadius(timeProvider.time())?problem.rho().outerValue():
                                                                                     problem.rho().innerValue());
      for(auto row=decltype(localRHS.size()){0};row!=localRHS.size();++row)
      {
        auto value(fdivValue*(velocityValue*phi[row])*rho);
        value*=weight;
        localRHS[row]+=value;
      }
    }
    localRHS.unbind();
  }
}

}
}

#endif // DUNE_FEM_CORRECTIONVELOCITYRHSNONDIVERGENCEFREE_HH
