#ifndef DUNE_FEM_COUPLEDOPERATORWRAPPER_HH
#define DUNE_FEM_COUPLEDOPERATORWRAPPER_HH

#include <dune/fem/operator/common/operator.hh>

namespace Dune
{
namespace Fem
{

template<typename VelocityOperatorImp,typename CurvatureVelocityOperatorImp,typename InterfaceOperatorImp,
         typename InterfaceInverseOperatorImp,typename VelocityCurvatureOperatorImp>
class CoupledOperatorWrapper:public Operator<typename VelocityOperatorImp::DomainFunctionType,
                                             typename VelocityOperatorImp::DomainFunctionType>
{
  public:
  typedef VelocityOperatorImp VelocityOperatorType;
  typedef CurvatureVelocityOperatorImp CurvatureVelocityOperatorType;
  typedef InterfaceOperatorImp InterfaceOperatorType;
  typedef InterfaceInverseOperatorImp InterfaceInverseOperatorType;
  typedef VelocityCurvatureOperatorImp VelocityCurvatureOperatorType;
  typedef CoupledOperatorWrapper<VelocityOperatorType,CurvatureVelocityOperatorType,InterfaceOperatorType,InterfaceInverseOperatorType,
                                 VelocityCurvatureOperatorType> ThisType;

  typedef typename VelocityOperatorType::DomainFunctionType DomainFunctionType;
  typedef DomainFunctionType RangeFunctionType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef DomainSpaceType RangeSpaceType;

  // constructor
  explicit CoupledOperatorWrapper(const VelocityOperatorType& velocityOp,const CurvatureVelocityOperatorType& curvatureVelocityOp,
                                  const InterfaceOperatorType& interfaceOp,InterfaceInverseOperatorType& interfaceInvOp,
                                  const VelocityCurvatureOperatorType& velocityCurvatureOp,double gamma):
    velocityop_(velocityOp),curvaturevelocityop_(curvatureVelocityOp),interfaceop_(interfaceOp),interfaceinvop_(interfaceInvOp),
    velocitycurvatureop_(velocityCurvatureOp),gamma_(gamma)
  {}

  CoupledOperatorWrapper(const ThisType& )=delete;

  const DomainSpaceType& domainSpace() const
  {
    return velocityop_.domainSpace();
  }

  const RangeSpaceType& rangeSpace() const
  {
    return velocityop_.rangeSpace();
  }

  void assemble() const
  {}

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    // w = velocityOperator*u
    velocityop_(u,w);

    if(gamma_!=0.0)
    {
      // rhs = [ velocityCurvatureOperator*u ; 0 ]
      typename InterfaceOperatorType::DomainFunctionType rhs("rhs",interfaceop_.domainSpace());
      velocitycurvatureop_(u,rhs.template subDiscreteFunction<0>());
      rhs.template subDiscreteFunction<1>().clear();

      // solution = interfaceInvOperator(rhs)
      typename InterfaceOperatorType::DomainFunctionType solution("solution",interfaceop_.domainSpace());
      interfaceinvop_.apply(rhs,solution);

      // velocity = curvatureVelocityOperator*solution_curvature
      DomainFunctionType velocity("velocity",velocityop_.domainSpace());
      curvaturevelocityop_(solution.template subDiscreteFunction<0>(),velocity);

      // w += gamma*velocity
      w.axpy(gamma_,velocity);
    }
  }

  private:
  const VelocityOperatorType& velocityop_;
  const CurvatureVelocityOperatorType& curvaturevelocityop_;
  const InterfaceOperatorType& interfaceop_;
  InterfaceInverseOperatorType& interfaceinvop_;
  const VelocityCurvatureOperatorType& velocitycurvatureop_;
  const double gamma_;
};

}
}
#endif // DUNE_FEM_COUPLEDOPERATORWRAPPER_HH
