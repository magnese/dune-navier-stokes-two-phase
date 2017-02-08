#ifndef DUNE_FEM_MESHSMOOTHING_HH
#define DUNE_FEM_MESHSMOOTHING_HH

#include <iostream>

#include <dune/common/timer.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/solver/umfpacksolver.hh>

#include "smoothingoperator.hh"
#include "problemssmoothing.hh"

namespace Dune
{
namespace Fem
{

template<typename FluidStateImp>
class MeshSmoothing
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef MeshSmoothing<FluidStateType> ThisType;

  // define discrete function and problem
  typedef typename FluidStateType::BulkDisplacementDiscreteFunctionType DiscreteFunctionType;
  typedef ParallelepipedGeometry<FluidStateType> SmoothingProblemType;

  // constructor
  explicit MeshSmoothing(FluidStateType& fluidState):
    fluidstate_(fluidState),problem_(fluidState),
    coeff_(Parameter::getValidValue<double>("CoeffSmoothing",1.0,[](auto val){return val>=0.0;})),isenabled_(coeff_>0.0?true:false)
  {}

  MeshSmoothing(const ThisType& )=delete;

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"Smoothing Coefficient = "<<coeff_<<(isenabled_?"":" (WARNING: smooth disabled!)")<<std::endl;
  }

  void enable()
  {
    isenabled_=true;
    if(coeff_<=0.0)
    {
      coeff_=1.0;
      std::cout<<"WARNING: setting smoothing coefficient = "<<coeff_<<std::endl;
    }
  }

  void disable()
  {
    isenabled_=false;
  }

  bool isEnabled() const
  {
    return isenabled_;
  }

  void computeBulkDisplacement()
  {
    // update fluid state
    fluidstate_.update();
    if(isenabled_)
      apply();
    else
    {
      fluidstate_.bulkDisplacement().clear();
      fluidstate_.meshManager().mapper().setInterfaceDFInBulkDF(fluidstate_.displacement(),fluidstate_.bulkDisplacement());
    }
  }

  private:
  FluidStateType& fluidstate_;
  SmoothingProblemType problem_;
  double coeff_;
  bool isenabled_;

  void apply()
  {
    // assemble operator
    Timer timerAssemble(false);
    timerAssemble.start();
    typedef SmoothingOperator<DiscreteFunctionType> SmoothingOperatorType;
    SmoothingOperatorType op(fluidstate_.bulkDisplacementSpace(),coeff_);
    op.assemble();
    timerAssemble.stop();
    // assemble RHS
    timerAssemble.start();
    DiscreteFunctionType rhs("smoothing RHS",fluidstate_.bulkDisplacementSpace());
    rhs.clear();
    fluidstate_.meshManager().mapper().setInterfaceDFInBulkDF(fluidstate_.displacement(),rhs);
    timerAssemble.stop();
    // impose BC
    timerAssemble.start();
    problem_.bc().apply(std::ignore,rhs,op);
    timerAssemble.stop();
    // impose given displacement for the interface
    constexpr std::size_t localBlockSize(FluidStateType::DisplacementDiscreteSpaceType::localBlockSize);
    const auto numBlocks(fluidstate_.displacement().blocks());
    for(auto i=decltype(numBlocks){0};i!=numBlocks;++i)
      for(auto l=decltype(localBlockSize){0};l!=localBlockSize;++l)
      {
        const auto row((fluidstate_.meshManager().mapper().vtxInterface2Bulk(i))*localBlockSize+l);
        op.systemMatrix().matrix().clearRow(row);
        op.systemMatrix().matrix().set(row,row,1.0);
      }
    // solve
    Timer timerSolve(false);
    timerSolve.start();
    UMFPACKOp<DiscreteFunctionType,SmoothingOperatorType> invOp(op);
    invOp(rhs,fluidstate_.bulkDisplacement());
    timerSolve.stop();
    // print timers
    std::cout<<"Assemble mesh smoothing operator time: "<<timerAssemble.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Solve mesh smoothing time: "<<timerSolve.elapsed()<<" seconds."<<std::endl;
  }
};

}
}

#endif // DUNE_FEM_MESHSMOOTHING_HH
