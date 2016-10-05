#ifndef DUNE_FEM_ASSEMBLEVELOCITYRHS_HH
#define DUNE_FEM_ASSEMBLEVELOCITYRHS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>

#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionType,typename ProblemType,typename TimeProviderType>
void assembleVelocityRHS(DiscreteFunctionType& rhs,const DiscreteFunctionType& oldVelocity,const ProblemType& problem,
                         const TimeProviderType& timeProvider)
{
  typedef typename DiscreteFunctionType::LocalFunctionType::RangeType LocalFunctionRangeType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
  const auto& space(rhs.space());
  std::vector<LocalFunctionRangeType> phi(space.blockMapper().maxNumDofs()*localBlockSize);

  // perform a grid walkthrough and assemble the RHS
  for(const auto& entity:space)
  {
    auto localRHS(rhs.localFunction(entity));
    const auto localOldVelocity(oldVelocity.localFunction(entity));
    const auto& baseSet(space.basisFunctionSet(entity));
    const auto rho(problem.rho(entity));

    CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space.order()+1);
    for(const auto& qp:quadrature)
    {
      const auto fValue(problem.velocityRHS()(entity.geometry().global(qp.position()),timeProvider.time(),entity));
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
