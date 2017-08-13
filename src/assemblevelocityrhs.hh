#ifndef DUNE_FEM_ASSEMBLEVELOCITYRHS_HH
#define DUNE_FEM_ASSEMBLEVELOCITYRHS_HH

#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionType,typename FluidStateType,typename ProblemType,typename TimeProviderType>
void assembleVelocityRHS(DiscreteFunctionType& rhs,const FluidStateType& fluidState,ProblemType& problem,
                         const TimeProviderType& timeProvider,const DiscreteFunctionType& oldVelocity)
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
  const auto& space(rhs.space());
  std::vector<typename DiscreteFunctionType::RangeType> phi(space.maxNumDofs());

  // perform a grid walkthrough and assemble the RHS
  LocalContribution<DiscreteFunctionType,Assembly::Add> localRHS(rhs);
  ConstLocalDiscreteFunction<DiscreteFunctionType> localOldVelocity(oldVelocity);
  #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
  typedef typename FluidStateType::PhysicalCoefficientDiscreteFunctionType PhysicalCoefficientDiscreteFunctionType;
  ConstLocalDiscreteFunction<PhysicalCoefficientDiscreteFunctionType> localOldRho(fluidState.rho());
  #endif
  for(const auto& entity:space)
  {
    localRHS.bind(entity);
    localOldVelocity.init(entity);
    const auto& baseSet(space.basisFunctionSet(entity));
    auto rho(problem.rho(entity));
    #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
    localOldRho.init(entity);
    typename PhysicalCoefficientDiscreteFunctionType::RangeType oldRho;
    localOldRho.evaluate(entity.geometry().center(),oldRho);
    rho=oldRho[0];
    #endif

    const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space.order()+1);
    for(const auto& qp:quadrature)
    {
      baseSet.evaluateAll(qp,phi);
      const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

      const auto localSize(localRHS.size());
      for(auto row=decltype(localSize){0};row!=localSize;++row)
      {
        typename DiscreteSpaceType::RangeFieldType value(0.0);
        for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
          for(auto kk=decltype(localSize){0};kk!=localSize;++kk)
            value+=localOldVelocity[kk]*phi[kk][k]*phi[row][k];
        value*=weight*(rho/timeProvider.deltaT());
        localRHS[row]+=value;
      }
    }
    localRHS.unbind();
  }
}

}
}

#endif // DUNE_FEM_ASSEMBLEVELOCITYRHS_HH
