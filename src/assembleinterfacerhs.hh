#ifndef DUNE_FEM_ASSEMBLEINTERFACERHS_HH
#define DUNE_FEM_ASSEMBLEINTERFACERHS_HH

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionType,typename OperatorType>
void assembleInterfaceRHS(DiscreteFunctionType& rhs,const OperatorType& op)
{
  DiscreteFunctionType temp("temp",rhs.space());
  temp.template subDiscreteFunction<0>().clear();
  temp.template subDiscreteFunction<1>().assign(rhs.space().grid().coordFunction().discreteFunction());

  op(temp,rhs);

  rhs.template subDiscreteFunction<0>().clear();
  rhs.template subDiscreteFunction<1>()*=-1.0;
}

}
}

#endif // DUNE_FEM_ASSEMBLEINTERFACERHS_HH
