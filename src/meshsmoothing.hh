#ifndef DUNE_FEM_MESHSMOOTHING_HH
#define DUNE_FEM_MESHSMOOTHING_HH

#include <iostream>
#include <memory>

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

template<class CoupledMeshManagerImp>
class MeshSmoothing
{
  public:
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef MeshSmoothing<CoupledMeshManagerType> ThisType;

  // define grid, function, space and grid part
  typedef typename CoupledMeshManagerType::BulkGridType GridType;
  typedef LeafGridPart<GridType> GridPartType;
  typedef FunctionSpace<double,double,GridType::dimensionworld,GridType::dimensionworld> ContinuosSpaceType;
  typedef LagrangeDiscreteFunctionSpace<ContinuosSpaceType,GridPartType,1> DiscreteSpaceType;
  typedef ISTLBlockVectorDiscreteFunction<DiscreteSpaceType> DiscreteFunctionType;

  // define problem
  typedef ParallelepipedGeometry<DiscreteSpaceType,CoupledMeshManagerType> SmoothingProblemType;

  // constructor
  explicit MeshSmoothing(CoupledMeshManagerType& meshManager):
    meshmanager_(meshManager),problem_(meshmanager_),coeff_(Parameter::getValue<double>("CoeffSmoothing",1.0)),sequence_(0)
  {}

  MeshSmoothing(const ThisType& )=delete;

  inline void printInfo(std::ostream& s=std::cout) const
  {
    s<<"coeff_smooth = "<<coeff_<<(coeff_>0.0?"":" (WARNING: smooth disabled!)")<<std::endl;
  }

  template<class InterfaceFunction>
  const DiscreteFunctionType& operator()(const InterfaceFunction& interfaceDisplacement)
  {
    update();
    if(coeff_>0.0)
      smoothing(interfaceDisplacement);
    else
    {
      displacement_->clear();
      meshmanager_.mapper().addInterfaceDF2BulkDF(interfaceDisplacement,*displacement_);
    }
    return *displacement_;
  }

  private:
  CoupledMeshManagerType& meshmanager_;
  std::unique_ptr<GridPartType> gridpart_;
  std::unique_ptr<DiscreteSpaceType> space_;
  std::unique_ptr<DiscreteFunctionType> displacement_;
  SmoothingProblemType problem_;
  const double coeff_;
  unsigned int sequence_;

  void update()
  {
    // check if the mesh is changed
    if(sequence_!=meshmanager_.sequence())
    {
      // reset pointers to avoid danglig references
      displacement_.reset();
      space_.reset();
      gridpart_.reset();
      // create grid part, space and displacement
      gridpart_=std::unique_ptr<GridPartType>(new GridPartType(meshmanager_.bulkGrid()));
      space_=std::unique_ptr<DiscreteSpaceType>(new DiscreteSpaceType(*gridpart_));
      displacement_=std::unique_ptr<DiscreteFunctionType>(new DiscreteFunctionType("bulk displacement",*space_));
      // update sequence number
      sequence_=meshmanager_.sequence();
    }
  }

  template<class InterfaceFunction>
  void smoothing(const InterfaceFunction& interfaceDisplacement)
  {
    // create timers
    Timer timerAssemble(false);
    Timer timerSolve(false);
    // assemble operator and impose BC
    timerAssemble.start();
    typedef SparseRowLinearOperator<DiscreteFunctionType,DiscreteFunctionType> SmoothingLinearOperatorType;
    typedef SmoothingOperator<SmoothingLinearOperatorType> SmoothigOperatorType;
    SmoothigOperatorType op(*space_,coeff_);
    op.assemble();
    problem_.bc().applyToOperator(op);
    timerAssemble.stop();
    // assemble RHS and impose BC
    timerAssemble.start();
    typedef SmoothingRHS<DiscreteFunctionType> RHSType;
    RHSType RHS(*space_);
    RHS.assemble(interfaceDisplacement,meshmanager_.mapper());
    problem_.bc().applyToRHS(RHS.rhs());
    timerAssemble.stop();
    // impose zero displacement for the interface
    const auto localBlockSize(InterfaceFunction::DiscreteFunctionSpaceType::localBlockSize);
    const auto numBlocks(interfaceDisplacement.blocks());
    for(auto i=0;i!=numBlocks;++i)
      for(auto l=0;l!=localBlockSize;++l)
      {
        const auto row((meshmanager_.mapper().vtxInterface2Bulk(i))*localBlockSize+l);
        op.systemMatrix().matrix().clearRow(row);
        op.systemMatrix().matrix().set(row,row,1.0);
      }
    // solve
    timerSolve.start();
    typedef UMFPACKOp<DiscreteFunctionType,SmoothigOperatorType> SmoothingInverseOperatorType;
    SmoothingInverseOperatorType invOp(op);
    invOp(RHS.rhs(),*displacement_);
    timerSolve.stop();
    // print timers
    std::cout<<"Assemble mesh smoothing operator time: "<<timerAssemble.elapsed()<<" seconds."<<std::endl;
    std::cout<<"Solve mesh smoothing time: "<<timerSolve.elapsed()<<" seconds."<<std::endl;
  }
};

}
}

#endif // DUNE_FEM_MESHSMOOTHING_HH
