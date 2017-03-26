#ifndef DUNE_FEM_CURVATUREVELOCITYOPERATOR_HH
#define DUNE_FEM_CURVATUREVELOCITYOPERATOR_HH

#include <dune/fem/io/io.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <fstream>
#include <string>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DomainFunctionImp,typename RangeFunctionImp,typename CoupledMeshManagerImp,
         template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class CurvatureVelocityOperator:public Operator<DomainFunctionImp,RangeFunctionImp>
{
  public:
  typedef DomainFunctionImp DomainFunctionType;
  typedef RangeFunctionImp RangeFunctionType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef LinearOperatorImp<DomainFunctionImp,RangeFunctionType> LinearOperatorType;
  typedef DomainFunctionType CurvatureFunctionType;
  typedef RangeFunctionType VelocityFunctionType;
  typedef CurvatureVelocityOperator<DomainFunctionType,RangeFunctionType,CoupledMeshManagerType,LinearOperatorImp> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename CurvatureFunctionType::DiscreteFunctionSpaceType CurvatureSpaceType;
  typedef CurvatureSpaceType DomainSpaceType;
  typedef typename VelocityFunctionType::DiscreteFunctionSpaceType VelocitySpaceType;
  typedef VelocitySpaceType RangeSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  explicit CurvatureVelocityOperator(const CurvatureSpaceType& curvatureSpace,const VelocitySpaceType& velocitySpace,
                                     const CoupledMeshManagerType& meshManager):
    curvaturespace_(curvatureSpace),velocityspace_(velocitySpace),op_("curvature velocity operator",curvaturespace_,velocityspace_),
    meshmanager_(meshManager)
  {}

  CurvatureVelocityOperator(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  virtual void operator()(const CurvatureFunctionType& u,VelocityFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="curvature_velocity_matrix.dat",unsigned int offset=0) const
  {
    const std::string& path(Parameter::getValue<std::string>("fem.prefix","."));
    if(!directoryExists(path))
      createDirectory(path);
    std::ofstream ofs(path+"/"+filename);
    op_.matrix().print(ofs,offset);
  }

  const CurvatureSpaceType& domainSpace() const
  {
    return curvaturespace_;
  }

  const VelocitySpaceType& rangeSpace() const
  {
    return velocityspace_;
  }

  void assemble()
  {
    // allocate matrix
    Stencil<CurvatureSpaceType,VelocitySpaceType> stencil(curvaturespace_,velocityspace_);
    for(const auto& interfaceEntity:curvaturespace_)
      stencil.fill(interfaceEntity,meshmanager_.correspondingInnerBulkIntersection(interfaceEntity).inside());
    op_.reserve(stencil);

    // clear matrix
    op_.clear();

    // allocate basis vectors
    constexpr std::size_t curvatureLocalBlockSize(CurvatureSpaceType::localBlockSize);
    typedef typename CurvatureFunctionType::LocalFunctionType::RangeType CurvatureRangeType;
    std::vector<CurvatureRangeType> phiCurvature(curvaturespace_.blockMapper().maxNumDofs()*curvatureLocalBlockSize);
    constexpr std::size_t velocityLocalBlockSize(VelocitySpaceType::localBlockSize);
    typedef typename VelocityFunctionType::LocalFunctionType::RangeType VelocityRangeType;
    std::vector<VelocityRangeType> phiVelocity(velocityspace_.blockMapper().maxNumDofs()*velocityLocalBlockSize);

    // perform an interface walkthrough and assemble the global matrix
    for(const auto& interfaceEntity:curvaturespace_)
    {
      // extract the associated bulk intersection and bulk entity
      const auto intersection(meshmanager_.correspondingInnerBulkIntersection(interfaceEntity));
      const auto bulkEntity(intersection.inside());

      // compute normal to interface
      const auto normalVector(intersection.centerUnitOuterNormal());

      // extract local matrix
      auto localMatrix(op_.localMatrix(interfaceEntity,bulkEntity));

      // define basis functions
      const auto& curvatureBaseSet(localMatrix.domainBasisFunctionSet());
      const auto& velocityBaseSet(localMatrix.rangeBasisFunctionSet());

      // loop over quadrature nodes
      typedef CachingQuadrature<typename VelocitySpaceType::GridPartType,1> QuadratureType;
      const QuadratureType quadrature(velocityspace_.gridPart(),intersection,2*curvaturespace_.order()+1,QuadratureType::INSIDE);
      for(const auto& qp:quadrature)
      {
        curvatureBaseSet.evaluateAll(qp.localPosition(),phiCurvature);
        velocityBaseSet.evaluateAll(qp,phiVelocity);
        const auto weight(intersection.geometry().integrationElement(qp.localPosition())*qp.weight());

        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(auto i=decltype(rowLocalSize)(0);i!=rowLocalSize;++i)
        {
          for(auto j=decltype(columnLocalSize){0};j!=columnLocalSize;++j)
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
  LinearOperatorType op_;
  const CoupledMeshManagerType& meshmanager_;
};

}
}

#endif // DUNE_FEM_CURVATUREVELOCITYOPERATOR_HH
