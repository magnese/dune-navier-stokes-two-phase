#ifndef DUNE_FEM_BULKVELOCITYOPERATOR_HH
#define DUNE_FEM_BULKVELOCITYOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <string>
#include <fstream>
#include <vector>

#define USE_SYMMETRIC_FORMULATION 1

namespace Dune
{
namespace Fem
{

template<typename LinearOperatorImp,typename ProblemImp>
class BulkVelocityOperator:public Operator<typename LinearOperatorImp::DomainFunctionType,typename LinearOperatorImp::RangeFunctionType>
{
  public:
  typedef LinearOperatorImp LinearOperatorType;
  typedef ProblemImp ProblemType;
  typedef typename LinearOperatorType::DomainFunctionType DomainFunctionType;
  typedef typename LinearOperatorType::RangeFunctionType RangeFunctionType;
  typedef BulkVelocityOperator<LinearOperatorType,ProblemType> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef DomainFunctionType DiscreteFunctionType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef DomainSpaceType DiscreteSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  explicit BulkVelocityOperator(const DiscreteSpaceType& space,const ProblemType& problem,const DiscreteFunctionType& oldSolution):
    space_(space),op_("bulk velocity operator",space_,space_),problem_(problem),oldsolution_(oldSolution)
  {}

  BulkVelocityOperator(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="velocity_matrix.dat") const
  {
    std::ofstream ofs(filename);
    op_.matrix().print(ofs);
  }

  const DomainSpaceType& domainSpace() const
  {
    return space_;
  }

  const RangeSpaceType& rangeSpace() const
  {
    return space_;
  }

  template<typename TimeProviderType>
  void assemble(const TimeProviderType& timeProvider)
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
    for(const auto& entity:space_)
    {
      const auto localOldSolution(oldsolution_.localFunction(entity));
      auto localMatrix(op_.localMatrix(entity,entity));
      typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
      const auto& baseSet(localMatrix.domainBasisFunctionSet());
      const auto mu(problem_.mu(entity));
      const auto rho(problem_.rho(entity));

      CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto& qp:quadrature)
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
            #if USE_SYMMETRIC_FORMULATION
            for(auto k=0;k!=localBlockSize;++k)
            {
              for(auto kk=0;kk!=localBlockSize;++kk)
              {
                value+=gradphi[i][k][kk]*gradphi[j][k][kk]+gradphi[i][kk][k]*gradphi[j][kk][k];
                value+=gradphi[i][k][kk]*gradphi[j][kk][k]+gradphi[i][kk][k]*gradphi[j][k][kk];
              }
            }
            value*=0.5;
            #else
            for(auto k=0;k!=localBlockSize;++k)
              value+=gradphi[i][k]*gradphi[j][k];
            #endif
            value*=mu;
            // convective part
            RangeFieldType valueConvective(0.0);
            for(auto k=0;k!=localBlockSize;++k)
            {
              for(auto kk=0;kk!=localBlockSize;++kk)
              {
                for(auto l=0;l!=localSize;++l)
                {
                  valueConvective+=localOldSolution[l]*phi[l][kk]*gradphi[j][k][kk]*phi[i][k];
                  valueConvective-=localOldSolution[l]*phi[l][kk]*gradphi[i][k][kk]*phi[j][k];
                }
              }
            }
            valueConvective*=0.5*rho;
            // time dependent part
            RangeFieldType valueTime(0.0);
            for(auto k=0;k!=localBlockSize;++k)
              valueTime+=phi[j][k]*phi[i][k];
            valueTime*=(rho/timeProvider.deltaT());
            // add to the local matrix
            value+=valueConvective;
            value+=valueTime;
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
  const ProblemType& problem_;
  // solution at the previous time step
  const DiscreteFunctionType& oldsolution_;
};

}
}

#endif // DUNE_FEM_BULKVELOCITYOPERATOR_HH
