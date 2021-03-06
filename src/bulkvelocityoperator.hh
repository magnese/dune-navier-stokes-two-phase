#ifndef DUNE_FEM_BULKVELOCITYOPERATOR_HH
#define DUNE_FEM_BULKVELOCITYOPERATOR_HH

#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/io/io.hh>
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

template<typename DiscreteFunctionImp,typename ProblemImp,template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class BulkVelocityOperator:public Operator<DiscreteFunctionImp,DiscreteFunctionImp>
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef ProblemImp ProblemType;
  typedef DiscreteFunctionType DomainFunctionType;
  typedef DiscreteFunctionType RangeFunctionType;
  typedef LinearOperatorImp<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef BulkVelocityOperator<DiscreteFunctionType,ProblemType,LinearOperatorImp> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef DomainSpaceType DiscreteSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  explicit BulkVelocityOperator(const DiscreteSpaceType& space,const ProblemType& problem,const DiscreteFunctionType& oldVelocity):
    space_(space),op_("bulk velocity operator",space_,space_),problem_(problem),oldvelocity_(oldVelocity)
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

  void print(const std::string& filename="velocity_matrix.dat",unsigned int offset=0) const
  {
    const std::string& path(Parameter::getValue<std::string>("fem.prefix","."));
    if(!directoryExists(path))
      createDirectory(path);
    std::ofstream ofs(path+"/"+filename);
    op_.exportMatrix().print(ofs,offset);
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

    constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
    std::vector<typename DiscreteFunctionType::RangeType> phi(space_.maxNumDofs());
    std::vector<typename DiscreteFunctionType::JacobianRangeType> gradphi(space_.maxNumDofs());
    ConstLocalDiscreteFunction<DiscreteFunctionType> localOldVelocity(oldvelocity_);
    #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
    typedef typename ProblemType::FluidStateType::PhysicalCoefficientDiscreteFunctionType PhysicalCoefficientDiscreteFunctionType;
    ConstLocalDiscreteFunction<PhysicalCoefficientDiscreteFunctionType> localOldRho(problem_.fluidState().rho());
    #endif

    // perform a grid walkthrough and assemble the global matrix
    for(const auto& entity:space_)
    {
      localOldVelocity.init(entity);
      auto localMatrix(op_.localMatrix(entity,entity));
      typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
      const auto& baseSet(localMatrix.domainBasisFunctionSet());
      const auto mu(problem_.mu(entity));
      const auto rho(problem_.rho(entity));

      #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
      localOldRho.init(entity);
      typename PhysicalCoefficientDiscreteFunctionType::RangeType oldRho;
      localOldRho.evaluate(entity.geometry().center(),oldRho);
      #endif

      const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto& qp:quadrature)
      {
        baseSet.evaluateAll(qp,phi);
        baseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

        const auto localSize(localMatrix.rows());
        for(auto i=decltype(localSize){0};i!=localSize;++i)
          for(auto j=decltype(localSize){0};j!=localSize;++j)
          {
            // laplacian part
            RangeFieldType value(0.0);
            #if USE_SYMMETRIC_LAPLACIAN_TERM
            for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
              for(auto kk=decltype(localBlockSize){0};kk!=localBlockSize;++kk)
              {
                value+=gradphi[i][k][kk]*gradphi[j][k][kk]+gradphi[i][kk][k]*gradphi[j][kk][k];
                value+=gradphi[i][k][kk]*gradphi[j][kk][k]+gradphi[i][kk][k]*gradphi[j][k][kk];
              }
            value*=0.5;
            #else
            for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
              value+=gradphi[i][k]*gradphi[j][k];
            #endif
            value*=mu;
            // navier-stokes only
            if(!problem_.isDensityNull())
            {
              // convective part
              RangeFieldType valueConvective(0.0);
              #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
              for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
                for(auto kk=decltype(localBlockSize){0};kk!=localBlockSize;++kk)
                  for(auto l=decltype(localSize){0};l!=localSize;++l)
                  {
                    valueConvective+=localOldVelocity[l]*phi[l][kk]*gradphi[j][k][kk]*phi[i][k];
                    valueConvective-=localOldVelocity[l]*phi[l][kk]*gradphi[i][k][kk]*phi[j][k];
                  }
              valueConvective*=0.5;
              #else
              for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
                for(auto kk=decltype(localBlockSize){0};kk!=localBlockSize;++kk)
                  for(auto l=decltype(localSize){0};l!=localSize;++l)
                    valueConvective+=localOldVelocity[l]*phi[l][kk]*gradphi[j][k][kk]*phi[i][k];
              #endif
              valueConvective*=rho;
              value+=valueConvective;
              // time dependent part
              RangeFieldType valueTime(0.0);
              for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
                valueTime+=phi[j][k]*phi[i][k];
              #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
              valueTime*=((rho+oldRho[0])/(2*timeProvider.deltaT()));
              #else
              valueTime*=(rho/timeProvider.deltaT());
              #endif
              value+=valueTime;
            }
            // add to the local matrix
            value*=weight;
            localMatrix.add(i,j,value);
          }
      }
    }
  }

  template<typename TimeProviderType>
  void allocateAndAssembleTimeDerivative(const TimeProviderType& timeProvider)
  {
    DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    op_.reserve(stencil);
    op_.clear();

    // navier-stokes only
    if(!problem_.isDensityNull())
    {
      constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
      std::vector<typename DiscreteFunctionType::RangeType> phi(space_.maxNumDofs());
      #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
      typedef typename ProblemType::FluidStateType::PhysicalCoefficientDiscreteFunctionType PhysicalCoefficientDiscreteFunctionType;
      ConstLocalDiscreteFunction<PhysicalCoefficientDiscreteFunctionType> localOldRho(problem_.fluidState().rho());
      #endif

      // perform a grid walkthrough and assemble the time derivative
      for(const auto& entity:space_)
      {
        auto localMatrix(op_.localMatrix(entity,entity));
        typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
        const auto& baseSet(localMatrix.domainBasisFunctionSet());
        auto rho(problem_.rho(entity));

        #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
        localOldRho.init(entity);
        typename PhysicalCoefficientDiscreteFunctionType::RangeType oldRho;
        localOldRho.evaluate(entity.geometry().center(),oldRho);
        rho+=oldRho[0];
        rho*=0.5;
        #endif

        const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
        for(const auto& qp:quadrature)
        {
          baseSet.evaluateAll(qp,phi);
          const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

          const auto localSize(localMatrix.rows());
          for(auto i=decltype(localSize){0};i!=localSize;++i)
            for(auto j=decltype(localSize){0};j!=localSize;++j)
            {
              RangeFieldType value(0.0);
              for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
                value+=phi[j][k]*phi[i][k];
              value*=(rho*weight/timeProvider.deltaT());
              localMatrix.add(i,j,value);
            }
        }
      }
    }
  }

  template<typename TimeProviderType>
  void assembleRemainingTerms(const TimeProviderType& timeProvider)
  {
    constexpr std::size_t localBlockSize(DiscreteSpaceType::localBlockSize);
    std::vector<typename DiscreteFunctionType::RangeType> phi(space_.maxNumDofs());
    std::vector<typename DiscreteFunctionType::JacobianRangeType> gradphi(space_.maxNumDofs());
    ConstLocalDiscreteFunction<DiscreteFunctionType> localOldVelocity(oldvelocity_);

    // perform a grid walkthrough and assemble the remaining terms
    for(const auto& entity:space_)
    {
      localOldVelocity.init(entity);
      auto localMatrix(op_.localMatrix(entity,entity));
      typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
      const auto& baseSet(localMatrix.domainBasisFunctionSet());
      const auto mu(problem_.mu(entity));
      const auto rho(problem_.rho(entity));

      const CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto& qp:quadrature)
      {
        baseSet.evaluateAll(qp,phi);
        baseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

        const auto localSize(localMatrix.rows());
        for(auto i=decltype(localSize){0};i!=localSize;++i)
          for(auto j=decltype(localSize){0};j!=localSize;++j)
          {
            // laplacian part
            RangeFieldType value(0.0);
            #if USE_SYMMETRIC_LAPLACIAN_TERM
            for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
              for(auto kk=decltype(localBlockSize){0};kk!=localBlockSize;++kk)
              {
                value+=gradphi[i][k][kk]*gradphi[j][k][kk]+gradphi[i][kk][k]*gradphi[j][kk][k];
                value+=gradphi[i][k][kk]*gradphi[j][kk][k]+gradphi[i][kk][k]*gradphi[j][k][kk];
              }
            value*=0.5;
            #else
            for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
              value+=gradphi[i][k]*gradphi[j][k];
            #endif
            value*=mu;
            // navier-stokes only
            if(!problem_.isDensityNull())
            {
              // convective part
              RangeFieldType valueConvective(0.0);
              #if USE_ANTISYMMETRIC_CONVECTIVE_TERM
              for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
                for(auto kk=decltype(localBlockSize){0};kk!=localBlockSize;++kk)
                  for(auto l=decltype(localSize){0};l!=localSize;++l)
                  {
                    valueConvective+=localOldVelocity[l]*phi[l][kk]*gradphi[j][k][kk]*phi[i][k];
                    valueConvective-=localOldVelocity[l]*phi[l][kk]*gradphi[i][k][kk]*phi[j][k];
                  }
              valueConvective*=0.5;
              #else
              for(auto k=decltype(localBlockSize){0};k!=localBlockSize;++k)
                for(auto kk=decltype(localBlockSize){0};kk!=localBlockSize;++kk)
                  for(auto l=decltype(localSize){0};l!=localSize;++l)
                    valueConvective+=localOldVelocity[l]*phi[l][kk]*gradphi[j][k][kk]*phi[i][k];
              #endif
              valueConvective*=rho;
              value+=valueConvective;
            }
            // add to the local matrix
            value*=weight;
            localMatrix.add(i,j,value);
          }
      }
    }
  }

  private:
  const DiscreteSpaceType& space_;
  LinearOperatorType op_;
  const ProblemType& problem_;
  // velocity at the previous step
  const DiscreteFunctionType& oldvelocity_;
};

}
}

#endif // DUNE_FEM_BULKVELOCITYOPERATOR_HH
