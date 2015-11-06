#ifndef DUNE_FEM_BULKVELOCITYRHS_HH
#define DUNE_FEM_BULKVELOCITYRHS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>

#include <cmath>
#include <fstream>
#include <vector>
#include <string>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionImp,typename ProblemImp>
class BulkVelocityRHS
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef ProblemImp ProblemType;
  typedef BulkVelocityRHS<DiscreteFunctionType,ProblemType> ThisType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;

  explicit BulkVelocityRHS(const DiscreteSpaceType& space,const ProblemType& problem,const DiscreteFunctionType& oldSolution):
    space_(space),rhs_("velocity RHS",space_),problem_(problem),oldsolution_(oldSolution)
  {}

  BulkVelocityRHS(const ThisType& )=delete;

  inline DiscreteFunctionType& rhs() const
  {
    return rhs_;
  }

  // dump rhs vector into file
  inline void print(const std::string& filename="velocity_rhs.dat") const
  {
    std::ofstream ofs(filename);
    rhs_.print(ofs);
    ofs.close();
  }

  inline double norm() const
  {
    return sqrt(rhs_.scalarProductDofs(rhs_));
  }

  template<typename FunctionType,typename TimeProviderType>
  void assemble(const FunctionType& function,const TimeProviderType& timeProvider) const
  {
    rhs_.clear();

    typedef typename DiscreteFunctionType::LocalFunctionType::RangeType LocalFunctionRangeType;
    const auto localBlockSize(DiscreteSpaceType::localBlockSize);
    std::vector<LocalFunctionRangeType> phi(space_.blockMapper().maxNumDofs()*localBlockSize);

    // perform a grid walkthrough and assemble the global matrix
    for(const auto entity:space_)
    {
      auto localRHS(rhs_.localFunction(entity));
      const auto localOldSolution(oldsolution_.localFunction(entity));
      const auto& baseSet(space_.basisFunctionSet(entity));
      const auto rho(problem_.rho(entity));

      CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> pointSet(entity,2*space_.order()+1);
      for(auto pt=0;pt!=pointSet.nop();++pt)
      {
        const auto fValue(function(entity.geometry().global(pointSet.point(pt)),timeProvider.time(),entity));
        baseSet.evaluateAll(pointSet.point(pt),phi);
        const auto weight(entity.geometry().integrationElement(pointSet.point(pt))*pointSet.weight(pt));

        const auto numLocalBlocks(localRHS.numScalarDofs());
        const auto localSize(numLocalBlocks*localBlockSize);
        std::size_t row(0);
        for(auto localIdx=0;localIdx!=numLocalBlocks;++localIdx)
        {
          for(auto l=0;l!=localBlockSize;++l,++row)
          {
            auto value(fValue*phi[row]);

            typename DiscreteSpaceType::RangeFieldType temp(0.0);
            for(auto k=0;k!=localBlockSize;++k)
              for(auto kk=0;kk!=localSize;++kk)
                temp+=localOldSolution[kk]*phi[kk][k]*phi[row][k];
            temp*=(rho/timeProvider.deltaT());

            value+=temp;
            value*=weight;
            localRHS[row]+=value;
          }
        }
      }
    }
  }

  private:
  const DiscreteSpaceType& space_;
  mutable DiscreteFunctionType rhs_;
  const ProblemType& problem_;
  // solution at the previous time step
  const DiscreteFunctionType& oldsolution_;
};

}
}

#endif // DUNE_FEM_BULKVELOCITYRHS_HH
