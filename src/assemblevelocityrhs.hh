#ifndef DUNE_FEM_ASSEMBLEVELOCITYRHS_HH
#define DUNE_FEM_ASSEMBLEVELOCITYRHS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>

#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionType,typename FluidStateType,typename ProblemType,typename TimeProviderType>
void assembleVelocityRHS(DiscreteFunctionType& rhs,const FluidStateType& fluidState,ProblemType& problem,
                         const TimeProviderType& timeProvider)
{
  typedef typename DiscreteFunctionType::LocalFunctionType::RangeType LocalFunctionRangeType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
  const auto& space(rhs.space());
  std::vector<LocalFunctionRangeType> phi(space.blockMapper().maxNumDofs()*localBlockSize);
  problem.velocityRHS().initialize(timeProvider.time(),timeProvider.time());

  // perform a grid walkthrough and assemble the RHS
  for(const auto& entity:space)
  {
    problem.velocityRHS().init(entity);
    auto localRHS(rhs.localFunction(entity));
    const auto localOldVelocity(fluidState.velocity().localFunction(entity));
    const auto& baseSet(space.basisFunctionSet(entity));
    auto rho(problem.rho(entity));
    #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
    const auto localOldRho(fluidState.rho().localFunction(entity));
    typename FluidStateType::PhysicalCoefficientDiscreteFunctionType::RangeType oldRho;
    localOldRho.evaluate(entity.geometry().center(),oldRho);
    rho+=oldRho[0];
    rho*=0.5;
    #endif

    const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space.order()+1);
    for(const auto& qp:quadrature)
    {
      typename FluidStateType::VelocityDiscreteFunctionType::RangeType fValue;
      problem.velocityRHS().evaluate(qp.position(),fValue);
      baseSet.evaluateAll(qp,phi);
      const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

      const auto localSize(localRHS.size());
      for(auto row=decltype(localSize){0};row!=localSize;++row)
      {
        auto value(fValue*phi[row]);
        if(!problem.isDensityNull())
        {
          typename DiscreteSpaceType::RangeFieldType temp(0.0);
          for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
            for(auto kk=decltype(localSize){0};kk!=localSize;++kk)
              temp+=localOldVelocity[kk]*phi[kk][k]*phi[row][k];
          temp*=(rho/timeProvider.deltaT());
          value+=temp;
        }
        value*=weight;
        localRHS[row]+=value;
      }
    }
  }
}

}
}

#endif // DUNE_FEM_ASSEMBLEVELOCITYRHS_HH
