#ifndef DUNE_FEM_BULKPRESSURERHS_HH
#define DUNE_FEM_BULKPRESSURERHS_HH

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <cmath>
#include <fstream>
#include <vector>
#include <string>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionImp>
class BulkPressureRHS
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef BulkPressureRHS<DiscreteFunctionType> ThisType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;

  explicit BulkPressureRHS(const DiscreteSpaceType& space):
    space_(space),rhs_("pressure RHS",space_)
  {}

  BulkPressureRHS(const ThisType& )=delete;

  DiscreteFunctionType& rhs() const
  {
    return rhs_;
  }

  // dump rhs vector into file
  void print(const std::string& filename="pressure_rhs.dat") const
  {
    std::ofstream ofs(filename);
    rhs_.print(ofs);
  }

  double norm() const
  {
    return sqrt(rhs_.scalarProductDofs(rhs_));
  }

  template<typename BC,typename TimeProviderType>
  void assemble(const BC& bc,const TimeProviderType& timeProvider) const
  {
    // clear RHS
    rhs_.clear();

    // create discrete function to interpolate BC
    typedef ISTLBlockVectorDiscreteFunction<typename BC::DomainSpaceType> BCDiscreteFunctionType;
    BCDiscreteFunctionType g("g",bc.domainSpace());

    // compute domain volume and \int_{\partial\Omega} I\vec g . \vec n
    typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
    const auto& gridPart(space_.gridPart());
    RangeFieldType vol(0.0);
    RangeFieldType integral(0.0);
    for(const auto& entity:space_)
    {
      vol+=std::abs(entity.geometry().volume());
      if(entity.hasBoundaryIntersections())
      {
        auto gLocal(g.localFunction(entity));
        for(const auto& intersection:intersections(static_cast<typename DiscreteFunctionType::GridPartType::GridViewType>(gridPart),entity))
          if(intersection.boundary())
          {
            bc.localInterpolateBoundaryFunction(timeProvider.time(),intersection,g);
            const auto normal(intersection.centerUnitOuterNormal());
            typedef CachingQuadrature<typename DiscreteFunctionType::GridPartType,1> QuadratureType;
            QuadratureType quadrature(gridPart,intersection,2*bc.domainSpace().order()+1,QuadratureType::INSIDE);
            for(const auto& qp:quadrature)
            {
              typename BCDiscreteFunctionType::RangeType gValue;
              gLocal.evaluate(qp.position(),gValue);
              const auto weight(intersection.geometry().integrationElement(qp.localPosition())*qp.weight());
              integral+=(normal*gValue)*weight;
            }
          }
      }
    }

    // if coeff is not 0 assemble RHS which is coeff*(1,\phi)
    const auto coeff(integral/vol);
    if(std::abs(coeff)>1.e-14)
    {
      typedef typename DiscreteFunctionType::LocalFunctionType::RangeType LocalFunctionRangeType;
      std::vector<LocalFunctionRangeType> phi(space_.blockMapper().maxNumDofs()*DiscreteSpaceType::localBlockSize);
      for(const auto& entity:space_)
      {
        auto localRHS(rhs_.localFunction(entity));
        const auto& baseSet(space_.basisFunctionSet(entity));
        CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
        for(const auto& qp:quadrature)
        {
          baseSet.evaluateAll(qp,phi);
          const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
          for(auto row=decltype(localRHS.size()){0};row!=localRHS.size();++row)
            localRHS[row]+=phi[row]*weight*coeff;
        }
      }
    }
  }

  private:
  const DiscreteSpaceType& space_;
  mutable DiscreteFunctionType rhs_;
};

}
}

#endif // DUNE_FEM_BULKPRESSURERHS_HH
