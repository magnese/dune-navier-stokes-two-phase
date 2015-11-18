#ifndef DUNE_FEM_SMOOTHINGOPERATOR_HH
#define DUNE_FEM_SMOOTHINGOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <string>
#include <fstream>

namespace Dune
{
namespace Fem
{

template<typename LinearOperatorImp>
class SmoothingOperator:public Operator<typename LinearOperatorImp::DomainFunctionType,typename LinearOperatorImp::RangeFunctionType>
{
  public:
  typedef LinearOperatorImp LinearOperatorType;
  typedef typename LinearOperatorType::DomainFunctionType DomainFunctionType;
  typedef typename LinearOperatorType::RangeFunctionType RangeFunctionType;
  typedef SmoothingOperator<LinearOperatorType> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef DomainFunctionType DiscreteFunctionType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef DomainSpaceType DiscreteSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  explicit SmoothingOperator(const DiscreteSpaceType& space,const double& coeff):
    space_(space),op_("smoothing operator",space_,space_),coeff_(coeff)
  {}

  SmoothingOperator(const ThisType& )=delete;

  inline const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  // apply the operator
  virtual inline void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  // dump system matrix into file
  inline void print(const std::string& filename="smoothing_matrix.dat") const
  {
    std::ofstream ofs(filename);
    op_.matrix().print(ofs,1);
    ofs.close();
  }

  inline const DomainSpaceType& domainSpace() const
  {
    return space_;
  }

  inline const RangeSpaceType& rangeSpace() const
  {
    return space_;
  }

  void assemble() const
  {
    DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    op_.reserve(stencil);
    op_.clear();

    const auto localBlockSize(DiscreteSpaceType::localBlockSize);
    typedef typename DiscreteFunctionType::LocalFunctionType::RangeType LocalFunctionRangeType;
    std::vector<LocalFunctionRangeType> phi(space_.blockMapper().maxNumDofs()*localBlockSize);
    typedef typename DiscreteFunctionType::LocalFunctionType::JacobianRangeType LocalFunctionJacobianRangeType;
    std::vector<LocalFunctionJacobianRangeType> gradphi(space_.blockMapper().maxNumDofs()*localBlockSize);

    // perform a grid walkthrough and assemble the global matrix
    for(const auto entity:space_)
    {
      auto localMatrix(op_.localMatrix(entity,entity));
      typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
      const auto& baseSet(localMatrix.domainBasisFunctionSet());

      CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto qp:quadrature)
      {
        // evaluate the jacobians of all basis functions
        baseSet.evaluateAll(qp,phi);
        baseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

        const auto localSize(localMatrix.rows());
        for(auto i=0;i!=localSize;++i)
        {
          for(auto j=0;j!=localSize;++j)
          {
            // laplacian part
            RangeFieldType value(0.0);
            for(auto k=0;k!=localBlockSize;++k)
            {
              for(auto kk=0;kk!=localBlockSize;++kk)
              {
                value+=gradphi[i][k][kk]*gradphi[j][k][kk]+gradphi[i][kk][k]*gradphi[j][kk][k];
                value+=gradphi[i][k][kk]*gradphi[j][kk][k]+gradphi[i][kk][k]*gradphi[j][k][kk];
              }
            }
            value*=0.25;

            // divergence part
            RangeFieldType valueDivergenceTest(0.0);
            for(auto k=0;k!=localBlockSize;++k)
              valueDivergenceTest+=gradphi[i][k][k];
            RangeFieldType valueDivergenceTrial(0.0);
            for(auto k=0;k!=localBlockSize;++k)
              valueDivergenceTrial+=gradphi[j][k][k];

            // add to the local matrix
            value+=coeff_*valueDivergenceTest*valueDivergenceTrial;
            value*=weight;
            localMatrix.add(i,j,value);
          }
        }

      }
    }

  }

  private:
  const DiscreteSpaceType& space_;
  mutable LinearOperatorType op_;
  const double& coeff_;
};

}
}

#endif // DUNE_FEM_SMOOTHINGOPERATOR_HH
