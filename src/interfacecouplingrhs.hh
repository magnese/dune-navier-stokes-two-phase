#ifndef DUNE_FEM_INTERFACECOUPLINGRHS_HH
#define DUNE_FEM_INTERFACECOUPLINGRHS_HH

// deprecation warning
#warning ("WARNING : interfacecouplingrhs.hh is deprecated")

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <vector>
#include <cmath>

namespace Dune
{
namespace Fem
{

template<typename CurvatureDiscreteFunctionType,typename VelocityDiscreteFunctionType,typename CoupledMeshManagerType>
void addCouplingInterfaceRHS(CurvatureDiscreteFunctionType& rhs,const VelocityDiscreteFunctionType& velocitySolutiontm,
                             const CoupledMeshManagerType& meshManager)
{
  // define bulk grid, bulk grid part and interface grid
  const auto& bulkGrid(velocitySolutiontm.grid());
  const auto& bulkGridPart(velocitySolutiontm.gridPart());
  const auto& interfaceGrid(rhs.grid());

  // define space, local function and basis for velocity
  const auto& velocitySpace(velocitySolutiontm.space());
  typedef typename VelocityDiscreteFunctionType::LocalFunctionType::RangeType VelocityLocalFunctionRangeType;
  typedef typename VelocityDiscreteFunctionType::DiscreteFunctionSpaceType VelocitySpaceType;
  constexpr std::size_t velocityLocalBlockSize(VelocitySpaceType::localBlockSize);
  std::vector<VelocityLocalFunctionRangeType> phiVelocity(velocitySpace.blockMapper().maxNumDofs()*velocityLocalBlockSize);

  // define space, local function and basis for curvature
  const auto& curvatureSpace(rhs.space());
  typedef typename CurvatureDiscreteFunctionType::LocalFunctionType::RangeType CurvatureLocalFunctionRangeType;
  typedef typename CurvatureDiscreteFunctionType::DiscreteFunctionSpaceType CurvatureSpaceType;
  constexpr std::size_t curvatureLocalBlockSize(CurvatureSpaceType::localBlockSize);
  std::vector<CurvatureLocalFunctionRangeType> phiCurvature(curvatureSpace.blockMapper().maxNumDofs()*curvatureLocalBlockSize);

  // loop over interface entities
  for(const auto& interfaceEntity:curvatureSpace)
  {
    // extract the associated bulk intersection and bulk entity
    const auto intersection(meshManager.correspondingInnerBulkIntersection(interfaceEntity));
    const auto bulkEntity(intersection.inside());

    // create local functions for velocity and curvature
    auto localRHS(rhs.localFunction(interfaceEntity));
    auto localVelocity(velocitySolutiontm.localFunction(bulkEntity));

    // compute normal to interface entity
    const auto normalVector(intersection.centerUnitOuterNormal());

    // define basis functions for velocity and curvature
    const auto& velocityBaseSet(velocitySpace.basisFunctionSet(bulkEntity));
    const auto& curvatureBaseSet(curvatureSpace.basisFunctionSet(interfaceEntity));

    // loop over quadrature nodes
    typedef CachingQuadrature<typename VelocityDiscreteFunctionType::GridPartType,1> QuadratureType;
    QuadratureType quadrature(bulkGridPart,intersection,2*curvatureSpace.order()+1,QuadratureType::INSIDE);
    for(const auto& qp:quadrature)
    {
      curvatureBaseSet.evaluateAll(qp.localPosition(),phiCurvature);
      velocityBaseSet.evaluateAll(qp,phiVelocity);
      const auto weight(intersection.geometry().integrationElement(qp.localPosition())*qp.weight());

      const auto curvatureNumLocalBlocks(localRHS.numDofs()/CurvatureSpaceType::dimRange);
      const auto velocityNumLocalBlocks(localVelocity.numDofs()/VelocitySpaceType::dimRange);
      for(auto row=decltype(curvatureNumLocalBlocks){0};row!=(curvatureLocalBlockSize*curvatureNumLocalBlocks);++row)
      {
        auto value(phiCurvature[row][0]);

        typename CurvatureSpaceType::RangeFieldType temp(0.0);
        std::size_t k(0);
        for(auto localIdx=decltype(velocityNumLocalBlocks){0};localIdx!=velocityNumLocalBlocks;++localIdx)
        {
          for(auto l=decltype(velocityLocalBlockSize){0};l!=velocityLocalBlockSize;++l,++k)
            temp+=localVelocity[k]*phiVelocity[localIdx][l]*normalVector[l];
        }

        value*=(temp*weight);
        localRHS[row]+=value;
      }
    }
  }
}

}
}

#endif // DUNE_FEM_INTERFACECOUPLINGRHS_HH
