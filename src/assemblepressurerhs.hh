#ifndef DUNE_FEM_ASSEMBLEPRESSURERHS_HH
#define DUNE_FEM_ASSEMBLEPRESSURERHS_HH

#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <cmath>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionType,typename BoundaryConditionType,typename TimeProviderType>
void assemblePressureRHS(DiscreteFunctionType& rhs,BoundaryConditionType& bc,const TimeProviderType& timeProvider)
{
  // create necessary quantities in BC
  bc.update();

  // compute domain volume and \int_{\partial\Omega} I\vec g . \vec n
  typedef typename BoundaryConditionType::DiscreteSpaceType::RangeFieldType RangeFieldType;
  RangeFieldType vol(0.0);
  RangeFieldType integral(0.0);
  for(const auto& entity:bc.space())
  {
    vol+=std::abs(entity.geometry().volume());
    if(entity.hasBoundaryIntersections())
    {
      for(const auto& intersection:intersections(bc.space().gridPart(),entity))
        if(intersection.boundary())
        {
          const auto& gLocalDOFs(bc.localBoundaryDOFs(timeProvider.time(),intersection));
          typedef LocalFunction<typename BoundaryConditionType::DiscreteSpaceType::BasisFunctionSetType,
                                typename BoundaryConditionType::LocalBoundaryDOFsType> GLocalType;
          GLocalType gLocal(bc.space().basisFunctionSet(entity),gLocalDOFs);
          const auto normal(intersection.centerUnitOuterNormal());
          typedef CachingQuadrature<typename BoundaryConditionType::DiscreteSpaceType::GridPartType,1> QuadratureType;
          const QuadratureType quadrature(bc.space().gridPart(),intersection,2*bc.space().order()+1,QuadratureType::INSIDE);
          for(const auto& qp:quadrature)
          {
            typename BoundaryConditionType::RangeType gValue;
            gLocal.evaluate(qp.position(),gValue);
            const auto weight(intersection.geometry().integrationElement(qp.localPosition())*qp.weight());
            integral+=(normal*gValue)*weight;
          }
        }
    }
  }

  // if coeff is not 0 assemble RHS which is coeff*(1,\phi)
  const double nullTolerance(Parameter::getValue<double>("NullTolerance",1.e-12));
  const auto coeff(-integral/vol);
  if(std::abs(coeff)>nullTolerance)
  {
    LocalContribution<DiscreteFunctionType,Assembly::Add> localRHS(rhs);
    std::vector<typename DiscreteFunctionType::RangeType> phi(rhs.space().maxNumDofs());
    for(const auto& entity:rhs.space())
    {
      localRHS.bind(entity);
      const auto& baseSet(rhs.space().basisFunctionSet(entity));
      const CachingQuadrature<typename DiscreteFunctionType::GridPartType,0> quadrature(entity,2*rhs.space().order()+1);
      for(const auto& qp:quadrature)
      {
        baseSet.evaluateAll(qp,phi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        for(auto row=decltype(localRHS.size()){0};row!=localRHS.size();++row)
          localRHS[row]+=phi[row]*weight*coeff;
      }
      localRHS.unbind();
    }
  }
}

}
}

#endif // DUNE_FEM_ASSEMBLEPRESSURERHS_HH
