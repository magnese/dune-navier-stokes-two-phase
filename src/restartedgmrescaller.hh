#ifndef DUNE_FEM_RESTARTED_GMRES_CALLER_HH
#define DUNE_FEM_RESTARTED_GMRES_CALLER_HH

namespace Dune
{
namespace Fem
{

#include <dune/fem/io/parameter.hh>
#include <dune/istl/solvers.hh>

template<typename DiscreteFunctionType>
class RestartedGMResCaller
{
  public:
  template<typename OperatorType,typename PreconditionerType>
  RestartedGMResCaller(OperatorType& op,PreconditionerType& preconditioner):
    invop_(op,preconditioner,Parameter::getValue<double>("SolverTolerance",1.e-12),Parameter::getValue<int>("SolverRestart",5),
      Parameter::getValue<int>("SolverMaxIterations",1000),Parameter::getValue<int>("SolverVerbosity",0))
  {}

  void operator()(DiscreteFunctionType& x,const DiscreteFunctionType& rhs)
  {
    // copy rhs to keep it unmodify
    auto rhsCopy(rhs);
    // solve system
    InverseOperatorResult info;
    invop_.apply(x,rhsCopy,info);
  }

  private:
  Dune::RestartedGMResSolver<DiscreteFunctionType> invop_;
};

}
}

#endif // DUNE_FEM_RESTARTED_GMRES_CALLER_HH
