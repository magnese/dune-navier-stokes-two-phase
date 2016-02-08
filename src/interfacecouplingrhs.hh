#ifndef DUNE_FEM_INTERFACECOUPLINGRHS_HH
#define DUNE_FEM_INTERFACECOUPLINGRHS_HH

// deprecation warning
#warning ("WARNING : interfacecouplingrhs.hh is deprecated")

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <vector>
#include <cmath>
#include "normal.hh"

namespace Dune
{
namespace Fem
{

template<typename CurvatureDiscreteFunctionType,typename VelocityDiscreteFunctionType,typename BulkInterfaceGridMapperType>
void addCouplingInterfaceRHS(CurvatureDiscreteFunctionType& rhs,const VelocityDiscreteFunctionType& velocitySolutiontm,
                             const BulkInterfaceGridMapperType& mapper)
{
  // define bulk grid, bulk grid part and interface grid
  const auto& bulkGrid(velocitySolutiontm.grid());
  const auto& bulkGridPart(velocitySolutiontm.gridPart());
  const auto& interfaceGrid(rhs.grid());

  // define space, local function and basis for velocity
  const auto& velocitySpace(velocitySolutiontm.space());
  typedef typename VelocityDiscreteFunctionType::LocalFunctionType::RangeType VelocityLocalFunctionRangeType;
  const auto velocityLocalBlockSize(VelocityDiscreteFunctionType::DiscreteFunctionSpaceType::localBlockSize);
  std::vector<VelocityLocalFunctionRangeType> phiVelocity(velocitySpace.blockMapper().maxNumDofs()*velocityLocalBlockSize);

  // define space, local function and basis for curvature
  const auto& curvatureSpace(rhs.space());
  typedef typename CurvatureDiscreteFunctionType::LocalFunctionType::RangeType CurvatureLocalFunctionRangeType;
  const auto curvatureLocalBlockSize(CurvatureDiscreteFunctionType::DiscreteFunctionSpaceType::localBlockSize);
  std::vector<CurvatureLocalFunctionRangeType> phiCurvature(curvatureSpace.blockMapper().maxNumDofs()*curvatureLocalBlockSize);

  // define normal functor and normal vector
  constexpr auto worlddim(VelocityDiscreteFunctionType::GridType::dimensionworld);
  typedef typename CurvatureDiscreteFunctionType::GridType::ctype ctype;
  typedef Normal<ctype,worlddim> NormalType;
  NormalType normal;
  typename NormalType::NormalVectorType normalVector;

  // loop over interface entities
  for(const auto& interfaceEntity:curvatureSpace)
  {
    // extract the corresponding bulk grid entity
    const auto idxInterface(interfaceGrid.leafIndexSet().index(interfaceEntity));
    const auto& bulkEntity(bulkGrid.entity(mapper.entitySeedInterface2Bulk(idxInterface)));
    const auto& faceLocalIdx(mapper.faceLocalIdxInterface2Bulk(idxInterface));

    // extract the associated intersection
    auto intersectionIt(bulkGridPart.ibegin(bulkEntity));
    while(intersectionIt->indexInInside()!=faceLocalIdx)
      ++intersectionIt;
    const auto intersection(*intersectionIt);

    // create local functions for velocity and curvature
    auto localRHS(rhs.localFunction(interfaceEntity));
    auto localVelocity(velocitySolutiontm.localFunction(bulkEntity));

    // compute normal to interface entity
    normal(interfaceEntity,normalVector,faceLocalIdx);

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

      const auto curvatureNumLocalBlocks(localRHS.numScalarDofs());
      const auto velocityNumLocalBlocks(localVelocity.numScalarDofs());
      for(auto row=0;row!=(curvatureLocalBlockSize*curvatureNumLocalBlocks);++row)
      {
        auto value(phiCurvature[row][0]);

        typename CurvatureDiscreteFunctionType::DiscreteFunctionSpaceType::RangeFieldType temp(0.0);
        std::size_t k(0);
        for(auto localIdx=0;localIdx!=velocityNumLocalBlocks;++localIdx)
        {
          for(auto l=0;l!=velocityLocalBlockSize;++l,++k)
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
