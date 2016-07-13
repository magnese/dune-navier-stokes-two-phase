#ifndef DUNE_FEM_VELOCITYCURVATUREOPERATOR_HH
#define DUNE_FEM_VELOCITYCURVATUREOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include "normal.hh"

#include <string>
#include <fstream>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DomainFunctionImp,typename RangeFunctionImp,typename BulkInterfaceGridMapperImp,
         template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class VelocityCurvatureOperator:public Operator<DomainFunctionImp,RangeFunctionImp>
{
  public:
  typedef DomainFunctionImp DomainFunctionType;
  typedef RangeFunctionImp RangeFunctionType;
  typedef BulkInterfaceGridMapperImp BulkInterfaceGridMapperType;
  typedef LinearOperatorImp<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef DomainFunctionType VelocityFunctionType;
  typedef RangeFunctionType CurvatureFunctionType;
  typedef VelocityCurvatureOperator<DomainFunctionType,RangeFunctionType,BulkInterfaceGridMapperType,LinearOperatorImp> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename VelocityFunctionType::DiscreteFunctionSpaceType VelocitySpaceType;
  typedef VelocitySpaceType DomainSpaceType;
  typedef typename CurvatureFunctionType::DiscreteFunctionSpaceType CurvatureSpaceType;
  typedef CurvatureSpaceType RangeSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;
  typedef typename CurvatureSpaceType::GridType InterfaceGridType;
  typedef typename CurvatureSpaceType::GridPartType InterfaceGridPartType;
  typedef typename VelocitySpaceType::GridType BulkGridType;
  typedef typename VelocitySpaceType::GridPartType BulkGridPartType;

  explicit VelocityCurvatureOperator(const VelocitySpaceType& velocitySpace,const CurvatureSpaceType& curvatureSpace,
                                     const BulkInterfaceGridMapperType& mapper):
    curvaturespace_(curvatureSpace),velocityspace_(velocitySpace),op_("velocity curvature operator",velocityspace_,curvaturespace_),
    mapper_(mapper),interfacegrid_(curvaturespace_.grid()),interfacegridpart_(curvaturespace_.gridPart()),bulkgrid_(velocityspace_.grid()),
    bulkgridpart_(velocityspace_.gridPart())
  {}

  VelocityCurvatureOperator(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  virtual void operator()(const VelocityFunctionType& u,CurvatureFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="velocity_curvature_matrix.dat") const
  {
    std::ofstream ofs(filename);
    op_.matrix().print(ofs);
  }

  const VelocitySpaceType& domainSpace() const
  {
    return velocityspace_;
  }

  const CurvatureSpaceType& rangeSpace() const
  {
    return curvaturespace_;
  }

  void assemble()
  {
    // allocate matrix
    Stencil<VelocitySpaceType,CurvatureSpaceType> stencil(velocityspace_,curvaturespace_);
    for(const auto& interfaceEntity:curvaturespace_)
    {
      // fill stencil
      const auto interfaceIdx(interfacegrid_.leafIndexSet().index(interfaceEntity));
      const auto& bulkEntity(bulkgrid_.entity(mapper_.entitySeedInterface2Bulk(interfaceIdx)));
      stencil.fill(bulkEntity,interfaceEntity);
    }
    op_.reserve(stencil);

    // clear matrix
    op_.clear();

    // define local functions and basis
    constexpr std::size_t velocityLocalBlockSize(VelocitySpaceType::localBlockSize);
    typedef typename VelocityFunctionType::LocalFunctionType::RangeType VelocityRangeType;
    std::vector<VelocityRangeType> phiVelocity(velocityspace_.blockMapper().maxNumDofs()*velocityLocalBlockSize);
    constexpr std::size_t curvatureLocalBlockSize(CurvatureSpaceType::localBlockSize);
    typedef typename CurvatureFunctionType::LocalFunctionType::RangeType CurvatureRangeType;
    std::vector<CurvatureRangeType> phiCurvature(curvaturespace_.blockMapper().maxNumDofs()*curvatureLocalBlockSize);

    // define normal functor and normal vector
    constexpr unsigned int worlddim(BulkGridType::dimensionworld);
    typedef Normal<typename InterfaceGridType::ctype,worlddim> NormalType;
    NormalType normal;
    typename NormalType::NormalVectorType normalVector;

    // perform an interface walkthrough and assemble the global matrix
    for(const auto& interfaceEntity:curvaturespace_)
    {
      // extract the corresponding bulk entity
      const auto interfaceIdx(interfacegrid_.leafIndexSet().index(interfaceEntity));
      const auto bulkEntity(bulkgrid_.entity(mapper_.entitySeedInterface2Bulk(interfaceIdx)));
      const auto& faceLocalIdx(mapper_.faceLocalIdxInterface2Bulk(interfaceIdx));

      // extract the associated bulk intersection
      auto intersectionIt(bulkgridpart_.ibegin(bulkEntity));
      while(static_cast<std::size_t>(intersectionIt->indexInInside())!=faceLocalIdx)
        ++intersectionIt;
      const auto intersection(*intersectionIt);

      // compute normal to interface
      normal(interfaceEntity,normalVector,faceLocalIdx);

      // extract local matrix
      auto localMatrix(op_.localMatrix(bulkEntity,interfaceEntity));

      // define basis functions
      const auto& velocityBaseSet(localMatrix.domainBasisFunctionSet());
      const auto& curvatureBaseSet(localMatrix.rangeBasisFunctionSet());

      // loop over quadrature nodes
      typedef CachingQuadrature<BulkGridPartType,1> QuadratureType;
      QuadratureType quadrature(bulkgridpart_,intersection,2*curvaturespace_.order()+1,QuadratureType::INSIDE);
      for(const auto& qp:quadrature)
      {
        velocityBaseSet.evaluateAll(qp,phiVelocity);
        curvatureBaseSet.evaluateAll(qp.localPosition(),phiCurvature);
        const auto weight(intersection.geometry().integrationElement(qp.localPosition())*qp.weight());

        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(auto i=decltype(rowLocalSize){0};i!=rowLocalSize;++i)
        {
          for(auto j=decltype(columnLocalSize){0};j!=columnLocalSize;++j)
          {
            const auto value((phiVelocity[j]*normalVector)*phiCurvature[i]*weight);
            localMatrix.add(i,j,value);
          }
        }

      }
    }
  }

  private:
  const CurvatureSpaceType& curvaturespace_;
  const VelocitySpaceType& velocityspace_;
  LinearOperatorType op_;
  const BulkInterfaceGridMapperType& mapper_;
  const InterfaceGridType& interfacegrid_;
  const InterfaceGridPartType& interfacegridpart_;
  const BulkGridType& bulkgrid_;
  const BulkGridPartType& bulkgridpart_;
};

}
}

#endif // DUNE_FEM_VELOCITYCURVATUREOPERATOR_HH
