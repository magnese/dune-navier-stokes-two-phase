#ifndef DUNE_FEM_PHYSICALCOEFFICIENT_HH
#define DUNE_FEM_PHYSICALCOEFFICIENT_HH

#include <iostream>
#include <string>

#include <dune/fem/function/common/localfunctionadapter.hh>
#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionSpaceImpl,typename FluidStateImp>
class PhysicalCoefficient
{
  public:
  typedef DiscreteFunctionSpaceImpl DiscreteFunctionSpaceType;
  typedef FluidStateImp FluidStateType;
  typedef PhysicalCoefficient<DiscreteFunctionSpaceType,FluidStateType> ThisType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  template<typename Validator>
  PhysicalCoefficient(const FluidStateType& fluidState,std::string&& parameterName,const Validator& validator):
    fluidstate_(fluidState),name_(parameterName),innercoeff_(Parameter::getValidValue<double>(name_+"Inner",1.0,validator)),
    outercoeff_(Parameter::getValidValue<double>(name_+"Outer",1.0,validator)),isnull_((innercoeff_==0.0&&outercoeff_==0.0)?true:false),
    delta_(outercoeff_-innercoeff_)
  {}

  PhysicalCoefficient(const ThisType& )=default;
  PhysicalCoefficient(ThisType&& )=default;
  ThisType& operator=(const ThisType& )=default;
  ThisType& operator=(ThisType&& )=default;

  RangeType operator()(const EntityType& entity) const
  {
    return fluidstate_.meshManager().bulkInnerIndicatorFunction().contains(entity)?innercoeff_:outercoeff_;
  }

  template<class PointType>
  void evaluate(const PointType& ,RangeType& ret) const
  {
    ret=(*this)(entity());
  }

  template<class PointType>
  void jacobian(const PointType& ,JacobianRangeType& ret) const
  {
    ret=JacobianRangeType(0.0);
  }

  template<class PointType>
  void hessian(const PointType& ,HessianRangeType& ret) const
  {
    ret=HessianRangeType(0.0);
  }

  void init(const EntityType& entity) const
  {
    entity_=&entity;
  }

  template<typename Arg1,typename Arg2>
  void initialize(Arg1&& ,Arg2&& ) const
  {}

  const EntityType& entity() const
  {
    assert(entity_);
    return *entity_;
  }

  bool isNull() const
  {
    return isnull_;
  }

  double delta() const
  {
    return delta_;
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<name_<<"Inner = "<<innercoeff_<<std::endl;
    s<<name_<<"Outer = "<<outercoeff_<<std::endl;
  }

  std::string name() const
  {
    return name_;
  }

  private:
  mutable EntityType const* entity_;
  const FluidStateType& fluidstate_;
  const std::string name_;
  const double innercoeff_;
  const double outercoeff_;
  const bool isnull_;
  const double delta_;
};

}
}

#endif // DUNE_FEM_PHYSICALCOEFFICIENT_HH
