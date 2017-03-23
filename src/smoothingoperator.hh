#ifndef DUNE_FEM_SMOOTHINGOPERATOR_HH
#define DUNE_FEM_SMOOTHINGOPERATOR_HH

#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <fstream>
#include <string>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionImp,template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class SmoothingOperator:public Operator<DiscreteFunctionImp,DiscreteFunctionImp>
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef DiscreteFunctionType DomainFunctionType;
  typedef DiscreteFunctionType RangeFunctionType;
  typedef LinearOperatorImp<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef SmoothingOperator<DiscreteFunctionType,LinearOperatorImp> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef DomainSpaceType DiscreteSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  explicit SmoothingOperator(const DiscreteSpaceType& space,double coeff):
    space_(space),op_("smoothing operator",space_,space_),coeff_(coeff)
  {}

  SmoothingOperator(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="smoothing_matrix.dat",unsigned int offset=0) const
  {
    std::ofstream ofs(Parameter::getValue<std::string>("fem.prefix",".")+"/"+filename);
    op_.matrix().print(ofs,offset);
  }

  const DomainSpaceType& domainSpace() const
  {
    return space_;
  }

  const RangeSpaceType& rangeSpace() const
  {
    return space_;
  }

  void assemble()
  {
    DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    op_.reserve(stencil);
    op_.clear();

    constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
    typedef typename DiscreteFunctionType::LocalFunctionType::JacobianRangeType LocalFunctionJacobianRangeType;
    std::vector<LocalFunctionJacobianRangeType> gradphi(space_.blockMapper().maxNumDofs()*localBlockSize);

    // perform a grid walkthrough and assemble the global matrix
    for(const auto& entity:space_)
    {
      auto localMatrix(op_.localMatrix(entity,entity));
      typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
      const auto& baseSet(localMatrix.domainBasisFunctionSet());

      const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto& qp:quadrature)
      {
        baseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

        const auto localSize(localMatrix.rows());
        for(auto i=decltype(localSize){0};i!=localSize;++i)
        {
          for(auto j=decltype(localSize){0};j!=localSize;++j)
          {
            // laplacian part
            RangeFieldType value(0.0);
            for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
            {
              for(auto kk=decltype(localBlockSize){0};kk!=localBlockSize;++kk)
              {
                value+=gradphi[i][k][kk]*gradphi[j][k][kk]+gradphi[i][kk][k]*gradphi[j][kk][k];
                value+=gradphi[i][k][kk]*gradphi[j][kk][k]+gradphi[i][kk][k]*gradphi[j][k][kk];
              }
            }
            value*=0.25;

            // divergence part
            RangeFieldType valueDivergenceTest(0.0);
            for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
              valueDivergenceTest+=gradphi[i][k][k];
            RangeFieldType valueDivergenceTrial(0.0);
            for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
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
  LinearOperatorType op_;
  const double coeff_;
};

}
}

#endif // DUNE_FEM_SMOOTHINGOPERATOR_HH
