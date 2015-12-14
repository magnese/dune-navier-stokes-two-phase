#ifndef DUNE_FEM_MESHSMOOTHING_HH
#define DUNE_FEM_MESHSMOOTHING_HH

#include <iostream>

#include <dune/common/timer.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/umfpacksolver.hh>

#include "smoothingoperator.hh"
#include "smoothingrhs.hh"
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

  // define space, discrete function and problem
  typedef typename FluidStateType::CoupledMeshManagerType CoupledMeshManagerType;
  typedef typename FluidStateType::BulkDisplacementDiscreteSpaceType DiscreteSpaceType;
  typedef typename FluidStateType::BulkDisplacementDiscreteFunctionType DiscreteFunctionType;
  typedef ParallelepipedGeometry<DiscreteSpaceType,CoupledMeshManagerType> SmoothingProblemType;

  // constructor
  explicit MeshSmoothing(FluidStateType& fluidState):
    fluidstate_(fluidState),problem_(fluidState.meshManager()),coeff_(Parameter::getValue<double>("CoeffSmoothing",1.0)),
    isenabled_(coeff_>0.0?true:false)
  {}

  MeshSmoothing(const ThisType& )=delete;

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"coeff_smooth = "<<coeff_<<(isenabled_?"":" (WARNING: smooth disabled!)")<<std::endl;
  }

  bool isEnabled() const
  {
    return isenabled_;
  }

  void apply()
  {
    // update fluid state
    fluidstate_.update();
    // create timers
    Timer timerAssemble(false);
    Timer timerSolve(false);
    // assemble operator and impose BC
    timerAssemble.start();
    typedef SparseRowLinearOperator<DiscreteFunctionType,DiscreteFunctionType> SmoothingLinearOperatorType;
    typedef SmoothingOperator<SmoothingLinearOperatorType> SmoothingOperatorType;
    SmoothingOperatorType op(fluidstate_.bulkDisplacementSpace(),coeff_);
    op.assemble();
    problem_.bc().applyToOperator(op);
    timerAssemble.stop();
    // assemble RHS and impose BC
    timerAssemble.start();
    typedef SmoothingRHS<DiscreteFunctionType> RHSType;
    RHSType RHS(fluidstate_.bulkDisplacementSpace());
    RHS.assemble(fluidstate_.displacement(),fluidstate_.meshManager().mapper());
    problem_.bc().applyToRHS(RHS.rhs());
    timerAssemble.stop();
    // impose zero displacement for the interface
    const auto localBlockSize(FluidStateType::DisplacementDiscreteSpaceType::localBlockSize);
    const auto numBlocks(fluidstate_.displacement().blocks());
    for(auto i=0;i!=numBlocks;++i)
      for(auto l=0;l!=localBlockSize;++l)
      {
        const auto row((fluidstate_.meshManager().mapper().vtxInterface2Bulk(i))*localBlockSize+l);
        op.systemMatrix().matrix().clearRow(row);
        op.systemMatrix().matrix().set(row,row,1.0);
      }
    // solve
    timerSolve.start();
    typedef UMFPACKOp<DiscreteFunctionType,SmoothingOperatorType> SmoothingInverseOperatorType;
    SmoothingInverseOperatorType invOp(op);
    invOp(RHS.rhs(),fluidstate_.bulkDisplacement());
    timerSolve.stop();
    // print timers
    std::cout<<"Assemble mesh smoothing operator time: "<<timerAssemble.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Solve mesh smoothing time: "<<timerSolve.elapsed()<<" seconds."<<std::endl;
  }

  private:
  FluidStateType& fluidstate_;
  SmoothingProblemType problem_;
  const double coeff_;
  const bool isenabled_;
};

}
}

#endif // DUNE_FEM_MESHSMOOTHING_HH
