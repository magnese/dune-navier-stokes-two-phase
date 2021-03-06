#ifndef DUNE_FEM_MESHSMOOTHING_HH
#define DUNE_FEM_MESHSMOOTHING_HH

#include <iostream>
#include <vector>

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
    s<<"Smoothing Coefficient = "<<coeff_<<(isenabled_?"":" (WARNING: smooth disabled!)")<<"\n";
    problem_.bc().printInfo(s);
  }

  void enable()
  {
    isenabled_=true;
    if(coeff_<=0.0)
    {
      coeff_=1.0;
      std::cout<<"WARNING: setting smoothing coefficient = "<<coeff_<<"\n";
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
      fluidstate_.meshManager().setInterfaceDFInBulkDF(fluidstate_.displacement(),fluidstate_.bulkDisplacement());
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
    const auto& space(fluidstate_.bulkDisplacementSpace());
    typedef SmoothingOperator<DiscreteFunctionType> SmoothingOperatorType;
    SmoothingOperatorType op(space,coeff_);
    op.assemble();
    timerAssemble.stop();
    // assemble RHS
    timerAssemble.start();
    DiscreteFunctionType rhs("smoothing RHS",space);
    rhs.clear();
    fluidstate_.meshManager().setInterfaceDFInBulkDF(fluidstate_.displacement(),rhs);
    timerAssemble.stop();
    // impose BC
    timerAssemble.start();
    problem_.bc().apply(std::ignore,rhs,op);
    timerAssemble.stop();
    // impose given displacement for the interface
    for(const auto& interfaceEntity:elements(fluidstate_.interfaceGridPart()))
    {
      const auto intersection(fluidstate_.meshManager().correspondingInnerBulkIntersection(interfaceEntity));
      const auto bulkEntity(intersection.inside());
      auto localMatrix(op.systemMatrix().localMatrix(bulkEntity,bulkEntity));
      const auto faceLocalIndex(intersection.indexInInside());
      std::vector<bool> globalBlockDofsFilter(space.blockMapper().numDofs(bulkEntity));
      space.blockMapper().onSubEntity(bulkEntity,faceLocalIndex,1,globalBlockDofsFilter);
      constexpr std::size_t localBlockSize(DiscreteFunctionType::DiscreteFunctionSpaceType::localBlockSize);
      const auto numLocalBlocks(space.basisFunctionSet(bulkEntity).size()/localBlockSize);
      std::size_t row(0);
      for(auto localIdx=decltype(numLocalBlocks){0};localIdx!=numLocalBlocks;++localIdx)
      {
        if(globalBlockDofsFilter[localIdx])
        {
          for(auto l=decltype(localBlockSize){0};l!=localBlockSize;++l,++row)
          {
            localMatrix.clearRow(row);
            localMatrix.set(row,row,1.0);
          }
        }
        else
          row+=localBlockSize;
      }
    }
    // solve
    Timer timerSolve(false);
    timerSolve.start();
    UMFPACKInverseOperator<DiscreteFunctionType,typename SmoothingOperatorType::MatrixType> invOp;
    invOp.bind(op.systemMatrix());
    invOp(rhs,fluidstate_.bulkDisplacement());
    timerSolve.stop();
    // print timers
    std::cout<<"Assemble mesh smoothing operator time: "<<timerAssemble.elapsed()<<" seconds.\n";
    std::cout<<"Solve mesh smoothing time: "<<timerSolve.elapsed()<<" seconds.\n";
  }
};

}
}

#endif // DUNE_FEM_MESHSMOOTHING_HH
