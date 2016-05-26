#ifndef DUNE_FEM_COUPLEDOPERATORWRAPPER_HH
#define DUNE_FEM_COUPLEDOPERATORWRAPPER_HH

#include <algorithm>

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
      // curvature = velocityCurvatureOperator*u
      const auto& curvatureSpace(curvaturevelocityop_.domainSpace());
      typename CurvatureVelocityOperatorType::CurvatureFunctionType curvature("curvature",curvatureSpace);
      velocitycurvatureop_(u,curvature);

      // copy curvature into interfaceRHS
      const auto& interfaceSpace(interfaceop_.domainSpace());
      typedef typename InterfaceOperatorType::DomainFunctionType InterfaceFunctionType;
      InterfaceFunctionType interfaceRHS("interface RHS",interfaceSpace);
      interfaceRHS.clear();
      std::copy(curvature.dbegin(),curvature.dend(),interfaceRHS.dbegin());

      // solution = interfaceInvOperator(interfaceRHS)
      InterfaceFunctionType solution("solution",interfaceSpace);
      interfaceinvop_.apply(interfaceRHS,solution);

      // copy solution into curvature
      std::copy_n(solution.dbegin(),curvature.size(),curvature.dbegin());

      // velocity = curvatureVelocityOperator*curvature
      const auto& velocitySpace(velocityop_.domainSpace());
      DomainFunctionType velocity("velocity",velocitySpace);
      curvaturevelocityop_(curvature,velocity);

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
