#ifndef DUNE_FEM_CURVATUREVELOCITYOPERATOR_HH
#define DUNE_FEM_CURVATUREVELOCITYOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>

#include "normal.hh"

#include <string>
#include <fstream>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename LinearOperatorImp,typename BulkInterfaceGridMapperImp>
class CurvatureVelocityOperator:public Operator<typename LinearOperatorImp::DomainFunctionType,
                                                typename LinearOperatorImp::RangeFunctionType>
{
  public:
  typedef LinearOperatorImp LinearOperatorType;
  typedef BulkInterfaceGridMapperImp BulkInterfaceGridMapperType;
  typedef typename LinearOperatorType::DomainFunctionType CurvatureFunctionType;
  typedef CurvatureFunctionType DomainFunctionType;
  typedef typename LinearOperatorType::RangeFunctionType VelocityFunctionType;
  typedef VelocityFunctionType RangeFunctionType;
  typedef CurvatureVelocityOperator<LinearOperatorType,BulkInterfaceGridMapperType> ThisType;
  typedef Operator<CurvatureFunctionType,VelocityFunctionType> BaseType;
  typedef typename CurvatureFunctionType::DiscreteFunctionSpaceType CurvatureSpaceType;
  typedef CurvatureSpaceType DomainSpaceType;
  typedef typename VelocityFunctionType::DiscreteFunctionSpaceType VelocitySpaceType;
  typedef VelocitySpaceType RangeSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;
  typedef typename CurvatureSpaceType::GridType InterfaceGridType;
  typedef typename CurvatureSpaceType::GridPartType InterfaceGridPartType;
  typedef typename VelocitySpaceType::GridType BulkGridType;
  typedef typename VelocitySpaceType::GridPartType BulkGridPartType;

  explicit CurvatureVelocityOperator(const CurvatureSpaceType& curvatureSpace,const VelocitySpaceType& velocitySpace,
                                     const BulkInterfaceGridMapperType& mapper):
    curvaturespace_(curvatureSpace),velocityspace_(velocitySpace),op_("curvature velocity operator",curvaturespace_,velocityspace_),
    mapper_(mapper),interfacegrid_(curvaturespace_.grid()),interfacegridpart_(curvaturespace_.gridPart()),bulkgrid_(velocityspace_.grid()),
    bulkgridpart_(velocityspace_.gridPart())
  {}

  CurvatureVelocityOperator(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  // apply the operator
  virtual void operator()(const CurvatureFunctionType& u,VelocityFunctionType& w) const
  {
    op_.apply(u,w);
  }

  // dump system matrix into file
  void print(const std::string& filename="curvature_velocity_matrix.dat") const
  {
    std::ofstream ofs(filename);
    const auto rows(op_.matrix().rows());
    auto count(decltype(rows){0});
    for(auto row=decltype(rows){0};row!=rows;++row)
    {
      while(count<(op_.matrix().numNonZeros()*(row+1)))
      {
        const auto entry(op_.matrix().realValue(count));
        const auto value(entry.first);
        const auto col(entry.second);
        if((std::abs(value)>1.e-13)&&(col>-1))
          ofs<<row+1<<" "<<col+1<<" "<<value<<std::endl;
        ++count;
      }
    }
  }

  const CurvatureSpaceType& domainSpace() const
  {
    return curvaturespace_;
  }

  const VelocitySpaceType& rangeSpace() const
  {
    return velocityspace_;
  }

  void assemble() const
  {
    // allocate matrix
    Stencil<CurvatureSpaceType,VelocitySpaceType> stencil(curvaturespace_,velocityspace_);
    for(const auto interfaceEntity:curvaturespace_)
    {
      // fill stencil
      const auto interfaceIdx(interfacegrid_.leafIndexSet().index(interfaceEntity));
      const auto& bulkEntity(bulkgrid_.entity(mapper_.entitySeedInterface2Bulk(interfaceIdx)));
      stencil.fill(interfaceEntity,bulkEntity);
    }
    op_.reserve(stencil);

    // clear matrix
    op_.clear();

    // define local functions and basis
    const auto curvatureLocalBlockSize(CurvatureSpaceType::localBlockSize);
    typedef typename CurvatureFunctionType::LocalFunctionType::RangeType CurvatureRangeType;
    std::vector<CurvatureRangeType> phiCurvature(curvaturespace_.blockMapper().maxNumDofs()*curvatureLocalBlockSize);
    const auto velocityLocalBlockSize(VelocitySpaceType::localBlockSize);
    typedef typename VelocityFunctionType::LocalFunctionType::RangeType VelocityRangeType;
    std::vector<VelocityRangeType> phiVelocity(velocityspace_.blockMapper().maxNumDofs()*velocityLocalBlockSize);

    // define normal functor and normal vector
    constexpr auto worlddim(BulkGridType::dimensionworld);
    typedef Normal<typename InterfaceGridType::ctype,worlddim> NormalType;
    NormalType normal;
    typename NormalType::NormalVectorType normalVector;

    // perform an interface walkthrough and assemble the global matrix
    for(const auto interfaceEntity:curvaturespace_)
    {
      // extract the corresponding bulk entity
      const auto interfaceIdx(interfacegrid_.leafIndexSet().index(interfaceEntity));
      const auto& bulkEntity(bulkgrid_.entity(mapper_.entitySeedInterface2Bulk(interfaceIdx)));
      const auto& faceLocalIdx(mapper_.faceLocalIdxInterface2Bulk(interfaceIdx));

      // extract the associated bulk intersection
      auto intersectionIt(bulkgridpart_.ibegin(bulkEntity));
      while(intersectionIt->indexInInside()!=faceLocalIdx)
        ++intersectionIt;
      const auto intersection(*intersectionIt);

      // compute normal to interface
      normal(interfaceEntity,normalVector,faceLocalIdx);

      // extract local matrix
      auto localMatrix(op_.localMatrix(interfaceEntity,bulkEntity));

      // define basis functions
      const auto& curvatureBaseSet(localMatrix.domainBasisFunctionSet());
      const auto& velocityBaseSet(localMatrix.rangeBasisFunctionSet());

      // loop over quadrature nodes
      typedef CachingQuadrature<BulkGridPartType,1> QuadratureType;
      QuadratureType quadrature(bulkgridpart_,intersection,2*curvaturespace_.order()+1,QuadratureType::INSIDE);
      for(const auto qp:quadrature)
      {
        curvatureBaseSet.evaluateAll(qp.localPosition(),phiCurvature);
        velocityBaseSet.evaluateAll(qp,phiVelocity);
        const auto weight(intersection.geometry().integrationElement(qp.localPosition())*qp.weight());

        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(auto i=0;i!=rowLocalSize;++i)
        {
          for(auto j=0;j!=columnLocalSize;++j)
          {
            const auto value((phiVelocity[i]*normalVector)*phiCurvature[j]*weight);
            localMatrix.add(i,j,value);
          }
        }

      }
    }
  }

  private:
  const CurvatureSpaceType& curvaturespace_;
  const VelocitySpaceType& velocityspace_;
  mutable LinearOperatorType op_;
  const BulkInterfaceGridMapperType& mapper_;
  const InterfaceGridType& interfacegrid_;
  const InterfaceGridPartType& interfacegridpart_;
  const BulkGridType& bulkgrid_;
  const BulkGridPartType& bulkgridpart_;
};

}
}

#endif // DUNE_FEM_CURVATUREVELOCITYOPERATOR_HH
