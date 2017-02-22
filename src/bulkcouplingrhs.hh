#ifndef DUNE_FEM_BULKCOUPLINGRHS_HH
#define DUNE_FEM_BULKCOUPLINGRHS_HH

// deprecation warning
#warning ("WARNING : bulkcouplingrhs.hh is deprecated")

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <vector>
#include <cmath>
#include "normal.hh"

namespace Dune
{
namespace Fem
{

template<typename VelocityDiscreteFunctionType,typename CurvatureDiscreteFunctionType,typename BulkInterfaceGridMapperType>
void addCouplingBulkRHS(VelocityDiscreteFunctionType& rhs,double gamma,const CurvatureDiscreteFunctionType& curvatureSolutiontm,
                        const BulkInterfaceGridMapperType& mapper)
{
  if(gamma!=0.0)
  {
    // define bulk grid, bulk grid part and interface grid
    const auto& bulkGrid(rhs.grid());
    const auto& bulkGridPart(rhs.gridPart());
    const auto& interfaceGrid(curvatureSolutiontm.grid());

    // define space, local function and basis for velocity
    const auto& velocitySpace(rhs.space());
    typedef typename VelocityDiscreteFunctionType::LocalFunctionType::RangeType VelocityLocalFunctionRangeType;
    typedef typename VelocityDiscreteFunctionType::DiscreteFunctionSpaceType VelocitySpaceType;
    constexpr std::size_t velocityLocalBlockSize(VelocitySpaceType::localBlockSize);
    std::vector<VelocityLocalFunctionRangeType> phiVelocity(velocitySpace.blockMapper().maxNumDofs()*velocityLocalBlockSize);

    // define space, local function and basis for curvature
    const auto& curvatureSpace(curvatureSolutiontm.space());
    typedef typename CurvatureDiscreteFunctionType::LocalFunctionType::RangeType CurvatureLocalFunctionRangeType;
    typedef typename CurvatureDiscreteFunctionType::DiscreteFunctionSpaceType CurvatureSpaceType;
    constexpr std::size_t curvatureLocalBlockSize(CurvatureSpaceType::localBlockSize);
    std::vector<CurvatureLocalFunctionRangeType> phiCurvature(curvatureSpace.blockMapper().maxNumDofs()*curvatureLocalBlockSize);

    // loop over interface grid entities
    for(const auto& interfaceEntity:curvatureSpace)
    {
      // extract the corresponding bulk grid entity
      const auto idxInterface(interfaceGrid.leafIndexSet().index(interfaceEntity));
      const auto bulkEntity(bulkGrid.entity(mapper.entitySeedInterface2Bulk(idxInterface)));
      const auto& faceLocalIdx(mapper.faceLocalIdxInterface2Bulk(idxInterface));

      // extract the associated intersection
      auto intersectionIt(bulkGridPart.ibegin(bulkEntity));
      while(static_cast<std::size_t>(intersectionIt->indexInInside())!=faceLocalIdx)
        ++intersectionIt;
      const auto intersection(*intersectionIt);

      // create local functions for velocity and curvature
      auto localRHS(rhs.localFunction(bulkEntity));
      auto localCurvature(curvatureSolutiontm.localFunction(interfaceEntity));

      // compute normal to interface entity
      const auto normalVector(computeNormal(interfaceEntity,faceLocalIdx));

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

        const auto curvatureNumLocalBlocks(localCurvature.numDofs()/CurvatureSpaceType::dimRange);
        const auto velocityNumLocalBlocks(localRHS.numDofs()/VelocitySpaceType::dimRange);
        std::size_t row(0);
        for(auto localIdx=decltype(velocityNumLocalBlocks){0};localIdx!=velocityNumLocalBlocks;++localIdx)
        {
          for(auto l=decltype(velocityLocalBlockSize){0};l!=velocityLocalBlockSize;++l,++row)
          {
            auto value(phiVelocity[row]*normalVector);

            typename CurvatureSpaceType::RangeFieldType temp(0.0);
            for(auto k=decltype(curvatureNumLocalBlocks){0};k!=(curvatureLocalBlockSize*curvatureNumLocalBlocks);++k)
              temp+=localCurvature[k]*phiCurvature[k][0];

            value*=(temp*weight*gamma);
            localRHS[row]+=value;
          }
        }
      }
    }
  }
}

}
}

#endif // DUNE_FEM_BULKCOUPLINGRHS_HH
