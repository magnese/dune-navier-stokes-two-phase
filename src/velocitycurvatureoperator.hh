#ifndef DUNE_FEM_VELOCITYCURVATUREOPERATOR_HH
#define DUNE_FEM_VELOCITYCURVATUREOPERATOR_HH

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
class VelocityCurvatureOperator:public Operator<typename LinearOperatorImp::DomainFunctionType,
                                                typename LinearOperatorImp::RangeFunctionType>
{
  public:
  typedef LinearOperatorImp LinearOperatorType;
  typedef BulkInterfaceGridMapperImp BulkInterfaceGridMapperType;
  typedef typename LinearOperatorType::DomainFunctionType VelocityFunctionType;
  typedef VelocityFunctionType DomainFunctionType;
  typedef typename LinearOperatorType::RangeFunctionType CurvatureFunctionType;
  typedef CurvatureFunctionType RangeFunctionType;
  typedef VelocityCurvatureOperator<LinearOperatorType,BulkInterfaceGridMapperType> ThisType;
  typedef Operator<VelocityFunctionType,CurvatureFunctionType> BaseType;
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

  inline const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  // apply the operator
  virtual inline void operator()(const VelocityFunctionType& u,CurvatureFunctionType& w) const
  {
    op_.apply(u,w);
  }

  // dump system matrix into file
  void print(const std::string& filename="velocity_curvature_matrix.dat") const
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

  inline const VelocitySpaceType& domainSpace() const
  {
    return velocityspace_;
  }

  inline const CurvatureSpaceType& rangeSpace() const
  {
    return curvaturespace_;
  }

  void assemble() const
  {
    // allocate matrix
    Stencil<VelocitySpaceType,CurvatureSpaceType> stencil(velocityspace_,curvaturespace_);
    for(const auto interfaceEntity:curvaturespace_)
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
    const auto velocityLocalBlockSize(VelocitySpaceType::localBlockSize);
    typedef typename VelocityFunctionType::LocalFunctionType::RangeType VelocityRangeType;
    std::vector<VelocityRangeType> phiVelocity(velocityspace_.blockMapper().maxNumDofs()*velocityLocalBlockSize);
    const auto curvatureLocalBlockSize(CurvatureSpaceType::localBlockSize);
    typedef typename CurvatureFunctionType::LocalFunctionType::RangeType CurvatureRangeType;
    std::vector<CurvatureRangeType> phiCurvature(curvaturespace_.blockMapper().maxNumDofs()*curvatureLocalBlockSize);

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
      const auto bulkEntity(bulkgrid_.entity(mapper_.entitySeedInterface2Bulk(interfaceIdx)));
      const auto& faceLocalIdx(mapper_.faceLocalIdxInterface2Bulk(interfaceIdx));

      // extract the associated bulk intersection
      auto intersectionIt(bulkgridpart_.ibegin(bulkEntity));
      while(intersectionIt->indexInInside()!=faceLocalIdx)
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
      for(const auto qp:quadrature)
      {
        velocityBaseSet.evaluateAll(qp,phiVelocity);
        curvatureBaseSet.evaluateAll(qp.localPosition(),phiCurvature);
        const auto weight(intersection.geometry().integrationElement(qp.localPosition())*qp.weight());

        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(auto i=0;i!=rowLocalSize;++i)
        {
          for(auto j=0;j!=columnLocalSize;++j)
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
  mutable LinearOperatorType op_;
  const BulkInterfaceGridMapperType& mapper_;
  const InterfaceGridType& interfacegrid_;
  const InterfaceGridPartType& interfacegridpart_;
  const BulkGridType& bulkgrid_;
  const BulkGridPartType& bulkgridpart_;
};

}
}

#endif // DUNE_FEM_VELOCITYCURVATUREOPERATOR_HH
