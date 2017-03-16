#ifndef DUNE_FEM_PROBLEMS_HH
#define DUNE_FEM_PROBLEMS_HH

#include <cmath>
#include <functional>
#include <iostream>
#include <ostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/interpolate.hh>

#include "boundarycondition.hh"
#include "physicalcoefficient.hh"

namespace Dune
{
namespace Fem
{

// base problem
template<typename FluidStateImp,template<typename ,typename > typename... VelocityBCImp>
class BaseProblem
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef BaseProblem<FluidStateType,VelocityBCImp...> ThisType;

  typedef typename FluidStateType::BulkDiscreteSpaceType BulkDiscreteSpaceType;
  typedef typename FluidStateType::CoupledMeshManagerType CoupledMeshManagerType;

  typedef typename FluidStateType::VelocityDiscreteSpaceType VelocityDiscreteSpaceType;
  typedef typename VelocityDiscreteSpaceType::EntityType EntityType;
  typedef typename VelocityDiscreteSpaceType::DomainType VelocityDomainType;
  typedef typename VelocityDiscreteSpaceType::RangeType VelocityRangeType;
  typedef LocalAnalyticalFunctionBinder<VelocityDiscreteSpaceType> VelocityFunctionType;
  typedef LocalFunctionAdapter<VelocityFunctionType> AdaptedVelocityFunctionType;

  typedef typename FluidStateType::PressureDiscreteSpaceType PressureDiscreteSpaceType;
  typedef typename PressureDiscreteSpaceType::DomainType PressureDomainType;
  typedef typename PressureDiscreteSpaceType::RangeType PressureRangeType;
  typedef LocalAnalyticalFunctionBinder<PressureDiscreteSpaceType> PressureFunctionType;
  typedef LocalFunctionAdapter<PressureFunctionType> AdaptedPressureFunctionType;

  typedef std::tuple<VelocityBCImp<VelocityDiscreteSpaceType,CoupledMeshManagerType>...> VelocityBCsType;

  typedef typename FluidStateType::PhysicalCoefficientDiscreteSpaceType PhysicalCoefficientDiscreteSpaceType;
  typedef PhysicalCoefficient<PhysicalCoefficientDiscreteSpaceType,FluidStateType> PhysicalCoefficientType;
  typedef LocalFunctionAdapter<PhysicalCoefficientType> AdaptedPhysicalCoefficientType;

  BaseProblem(FluidStateType& fluidState,bool isTimeDependent,bool hasExactSolution,const std::string& name):
    fluidstate_(fluidState),velocitybcs_(VelocityBCImp<VelocityDiscreteSpaceType,CoupledMeshManagerType>(fluidstate_.meshManager())...),
    istimedependent_(isTimeDependent),hasexactsolution_(hasExactSolution),name_(name),mu_(fluidstate_,"Mu",[](auto val){return val>0.0;}),
    gamma_(Parameter::getValue<double>("Gamma",1.0)),rho_(fluidstate_,"Rho",[](auto val){return val>=0.0;})
  {}

  bool isTimeDependent() const
  {
    return istimedependent_;
  }
  bool hasExactSolution() const
  {
    return hasexactsolution_;
  }

  VelocityFunctionType& velocityIC()
  {
    return std::get<1>(velocity_);
  }
  const VelocityFunctionType& velocityIC() const
  {
    return std::get<1>(velocity_);
  }
  template<typename DF>
  void velocityIC(DF& df) const
  {
    if(istimedependent_)
    {
      AdaptedVelocityFunctionType velocityICAdapted("velocity IC adapted",velocityIC(),fluidstate_.bulkGridPart(),
        VelocityDiscreteSpaceType::polynomialOrder);
      velocityICAdapted.initialize(0.0,0.0);
      interpolate(velocityICAdapted,df);
    }
  }
  VelocityFunctionType& velocityRHS()
  {
    return std::get<2>(velocity_);
  }
  const VelocityFunctionType& velocityRHS() const
  {
    return std::get<2>(velocity_);
  }
  VelocityFunctionType& velocitySolution()
  {
    return std::get<0>(velocity_);
  }
  const VelocityFunctionType& velocitySolution() const
  {
    return std::get<0>(velocity_);
  }
  template<typename DF>
  void velocitySolution(DF& df,double t) const
  {
    if(hasexactsolution_)
    {
      AdaptedVelocityFunctionType velocitySolutionAdapted("velocity solution adapted",velocitySolution(),fluidstate_.bulkGridPart(),
        VelocityDiscreteSpaceType::polynomialOrder);
      velocitySolutionAdapted.initialize(t,t);
      interpolate(velocitySolutionAdapted,df);
    }
  }

  PressureFunctionType& pressureIC()
  {
    return std::get<1>(pressure_);
  }
  const PressureFunctionType& pressureIC() const
  {
    return std::get<1>(pressure_);
  }
  template<typename DF>
  void pressureIC(DF& df) const
  {
    if(istimedependent_)
    {
      AdaptedPressureFunctionType pressureICAdapted("pressure IC adapted",pressureIC(),fluidstate_.bulkGridPart(),
        PressureDiscreteSpaceType::polynomialOrder);
      pressureICAdapted.initialize(0.0,0.0);
      interpolate(pressureICAdapted,df);
    }
  }
  PressureFunctionType& pressureSolution()
  {
    return std::get<0>(pressure_);
  }
  const PressureFunctionType& pressureSolution() const
  {
    return std::get<0>(pressure_);
  }
  template<typename DF>
  void pressureSolution(DF& df,double t) const
  {
    if(hasexactsolution_)
    {
      AdaptedPressureFunctionType pressureSolutionAdapted("pressure solution adapted",pressureSolution(),fluidstate_.bulkGridPart(),
        PressureDiscreteSpaceType::polynomialOrder);
      pressureSolutionAdapted.initialize(t,t);
      interpolate(pressureSolutionAdapted,df);
    }
  }

  double mu(const EntityType& entity) const
  {
    return mu_(entity);
  }
  const PhysicalCoefficientType& mu() const
  {
    return mu_;
  }
  template<typename DF>
  void mu(DF& df) const
  {
    AdaptedPhysicalCoefficientType muAdapted(mu_.name()+" adapted",mu_,fluidstate_.bulkGridPart(),
      PhysicalCoefficientDiscreteSpaceType::polynomialOrder);
    interpolate(muAdapted,df);
  }
  double gamma() const
  {
    return gamma_;
  }
  double rho(const EntityType& entity) const
  {
    return rho_(entity);
  }
  const PhysicalCoefficientType& rho() const
  {
    return rho_;
  }
  template<typename DF>
  void rho(DF& df) const
  {
    AdaptedPhysicalCoefficientType rhoAdapted(rho_.name()+" adapted",rho_,fluidstate_.bulkGridPart(),
      PhysicalCoefficientDiscreteSpaceType::polynomialOrder);
    interpolate(rhoAdapted,df);
  }
  bool isDensityNull() const
  {
    return rho_.isNull();
  }

  auto& velocityBC()
  {
    return std::get<0>(velocitybcs_);
  }

  template<typename... Args>
  void applyBC(Args&... args)
  {
    Hybrid::forEach(std::make_index_sequence<std::tuple_size<VelocityBCsType>::value>{},
      [&](auto i){std::get<i>(velocitybcs_).apply(args...);});
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"Problem : "<<name_<<"\n";
    mu_.printInfo(s);
    s<<"Gamma = "<<gamma_<<"\n";
    rho_.printInfo(s);
    Hybrid::forEach(std::make_index_sequence<std::tuple_size<VelocityBCsType>::value>{},
      [&](auto i){std::get<i>(velocitybcs_).printInfo(s);});

  }

  protected:
  FluidStateType& fluidstate_;
  VelocityBCsType velocitybcs_;
  const bool istimedependent_;
  const bool hasexactsolution_;
  const std::string name_;
  const PhysicalCoefficientType mu_;
  const double gamma_;
  const PhysicalCoefficientType rho_;
  std::tuple<VelocityFunctionType,VelocityFunctionType,VelocityFunctionType> velocity_;
  std::tuple<PressureFunctionType,PressureFunctionType> pressure_;
  static constexpr unsigned int worlddim=CoupledMeshManagerType::BulkGridType::dimensionworld;

  template<std::size_t N>
  auto& getVelocityBC()
  {
    return std::get<N>(velocitybcs_);
  }
};

// Stokes test 1
template<typename FluidStateImp>
class StokesTest1Problem:public BaseProblem<FluidStateImp,DirichletCondition>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef StokesTest1Problem<FluidStateType> ThisType;
  typedef BaseProblem<FluidStateType,DirichletCondition> BaseType;

  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::VelocityDomainType VelocityDomainType;
  typedef typename BaseType::VelocityRangeType VelocityRangeType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::velocityRHS;
  using BaseType::mu;

  StokesTest1Problem(FluidStateType& fluidState):
    BaseType(fluidState,false,true,"Stokes test 1")
  {
    velocityRHS().function()=[&](const VelocityDomainType& x,double ,const EntityType& entity)
    {
      VelocityRangeType value(0.0);
      value[1]=-6.0*mu(entity)*x[0];
      return value;
    };

    velocitySolution().function()=[](const VelocityDomainType& x,double ,const EntityType& )
    {
      VelocityRangeType value(0.0);
      value[1]=std::pow(x[0],3)-x[0];
      return value;
    };

    velocityIC()=velocitySolution();

    velocityBC().addAllBoundaryIDs(velocitySolution());
  }
};

// Stokes test 2
template<typename FluidStateImp>
class StokesTest2Problem:public BaseProblem<FluidStateImp,DirichletCondition>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef StokesTest2Problem<FluidStateType> ThisType;
  typedef BaseProblem<FluidStateType,DirichletCondition> BaseType;

  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::VelocityDomainType VelocityDomainType;
  typedef typename BaseType::VelocityRangeType VelocityRangeType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::velocityRHS;
  using BaseType::mu;

  StokesTest2Problem(FluidStateType& fluidState):
    BaseType(fluidState,false,true,"Stokes test 2")
  {
    velocityRHS().function()=[&](const VelocityDomainType& x,double ,const EntityType& entity)
    {
      const auto muValue(mu(entity));
      VelocityRangeType value(0.0);
      value[0]=-12.0*muValue*std::pow(x[1],2);
      value[1]=-6.0*muValue*x[0];
      return value;
    };

    velocitySolution().function()=[](const VelocityDomainType& x,double ,const EntityType& )
    {
      VelocityRangeType value(0.0);
      value[0]=1.0+std::pow(x[1],4);
      value[1]=std::pow(x[0],3)-x[0];
      return value;
    };

    velocityIC()=velocitySolution();

    velocityBC().addAllBoundaryIDs(velocitySolution());
  }
};

// stationary bubble
template<typename FluidStateImp>
class StationaryBubbleProblem:public BaseProblem<FluidStateImp,DirichletCondition>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef StationaryBubbleProblem<FluidStateType> ThisType;
  typedef BaseProblem<FluidStateType,DirichletCondition> BaseType;

  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::PressureDomainType PressureDomainType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::pressureSolution;
  using BaseType::pressureIC;
  using BaseType::gamma;
  using BaseType::worlddim;
  using BaseType::fluidstate_;

  StationaryBubbleProblem(FluidStateType& fluidState):
    BaseType(fluidState,true,true,"stationary bubble"),r0_(0.5)
  {
    pressureSolution().function()=[&](const PressureDomainType& x,double ,const EntityType& )
    {
      const auto rt(exactRadius());
      const auto coeff(gamma()*static_cast<double>(worlddim-1)/rt);
      const auto indicatorValue(x.two_norm()<=rt?1.0:0.0);
      auto value(coeff*(indicatorValue-std::pow(4.0/3.0,worlddim-2)*M_PI*std::pow(rt,worlddim)/std::pow(2.0,worlddim)));
      return value;
    };

    pressureIC()=pressureSolution();

    velocityBC().addAllBoundaryIDs(velocitySolution());
  }

  template<typename... Args>
  double exactRadius(const Args&... ) const
  {
    return r0_;
  }

  private:
  const double r0_;
};

// expanding bubble
template<typename FluidStateImp>
class ExpandingBubbleProblem:public BaseProblem<FluidStateImp,DirichletCondition>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef ExpandingBubbleProblem<FluidStateType> ThisType;
  typedef BaseProblem<FluidStateType,DirichletCondition> BaseType;

  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::VelocityDomainType VelocityDomainType;
  typedef typename BaseType::PressureDomainType PressureDomainType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::pressureSolution;
  using BaseType::pressureIC;
  using BaseType::mu_;
  using BaseType::gamma;
  using BaseType::worlddim;
  using BaseType::fluidstate_;

  ExpandingBubbleProblem(FluidStateType& fluidState):
    BaseType(fluidState,true,true,"expanding bubble"),r0_(0.5),alpha_(0.15)
  {
    velocitySolution().function()=[&](const VelocityDomainType& x,double ,const EntityType& )
    {
      auto value(x);
      value*=alpha_;
      value/=std::pow(x.two_norm(),worlddim);
      return value;
    };

    pressureSolution().function()=[&](const PressureDomainType& x,double t,const EntityType& )
    {
      const auto rt(exactRadius(t));
      const auto coeff(static_cast<double>(worlddim-1)*(gamma()/rt+2.0*alpha_*mu_.delta()/std::pow(rt,worlddim)));
      const auto indicatorValue(x.two_norm()<=rt?1.0:0.0);
      auto value(coeff*(indicatorValue-(std::pow(4.0/3.0,worlddim-2)*M_PI*std::pow(rt,worlddim)-std::pow(2.0/3.0,worlddim))/
                                       (std::pow(2.0,worlddim)-std::pow(2.0/3.0,worlddim))));
      return value;
    };

    velocityIC()=velocitySolution();
    pressureIC()=pressureSolution();

    velocityBC().addAllBoundaryIDs(velocitySolution());
  }

  double exactRadius(double t) const
  {
    return std::pow(std::pow(r0_,worlddim)+alpha_*t*static_cast<double>(worlddim),1.0/static_cast<double>(worlddim));
  }

  protected:
  const double r0_;
  const double alpha_;
};

// shear flow
template<typename FluidStateImp>
class ShearFlowProblem:public BaseProblem<FluidStateImp,DirichletCondition>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef ShearFlowProblem<FluidStateType> ThisType;
  typedef BaseProblem<FluidStateType,DirichletCondition> BaseType;

  typedef typename BaseType::VelocityFunctionType VelocityFunctionType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::VelocityDomainType VelocityDomainType;
  typedef typename BaseType::VelocityRangeType VelocityRangeType;

  using BaseType::velocityBC;
  using BaseType::worlddim;

  ShearFlowProblem(FluidStateType& fluidState):
    BaseType(fluidState,true,false,"shear flow")
  {
    VelocityFunctionType f([&](const VelocityDomainType& x,double ,const EntityType& )
                           {
                             VelocityRangeType value(0.0);
                             value[0]=x[worlddim-1];
                             return value;
                           });

    velocityBC().addAllBoundaryIDs(f);
  }
};

// stationary Navier-Stokes
template<typename FluidStateImp>
class StationaryNavierStokesProblem:public BaseProblem<FluidStateImp,DirichletCondition>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef StationaryNavierStokesProblem<FluidStateType> ThisType;
  typedef BaseProblem<FluidStateType,DirichletCondition> BaseType;

  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::VelocityDomainType VelocityDomainType;
  typedef typename BaseType::VelocityRangeType VelocityRangeType;
  typedef typename BaseType::PressureDomainType PressureDomainType;
  typedef typename BaseType::PressureRangeType PressureRangeType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::velocityRHS;
  using BaseType::pressureSolution;
  using BaseType::pressureIC;
  using BaseType::mu;

  StationaryNavierStokesProblem(FluidStateType& fluidState):
    BaseType(fluidState,false,true,"stationary Navier-Stokes")
  {
    velocityRHS().function()=[&](const VelocityDomainType& x,double ,const EntityType& entity)
    {
      const auto muValue(mu(entity));
      VelocityRangeType value(0.0);
      value[0]=-2.0*std::pow(M_PI,2)*muValue*std::cos(M_PI*x[0])*std::sin(M_PI*x[1])+M_PI/2.0*std::sin(2.0*M_PI*x[0]);
      value[1]=2.0*std::pow(M_PI,2)*muValue*std::sin(M_PI*x[0])*std::cos(M_PI*x[1])+M_PI/2.0*std::sin(2.0*M_PI*x[1]);
      return value;
    };

    velocitySolution().function()=[](const VelocityDomainType& x,double ,const EntityType& )
    {
      VelocityRangeType value(0.0);
      value[0]=-std::cos(M_PI*x[0])*std::sin(M_PI*x[1]);
      value[1]=std::sin(M_PI*x[0])*std::cos(M_PI*x[1]);
      return value;
    };

    pressureSolution().function()=[](const PressureDomainType& x,double ,const EntityType& )
    {
      PressureRangeType value(0.0);
      value[0]=-0.25*(std::cos(2.0*M_PI*x[0])+std::cos(2.0*M_PI*x[1]));
      return value;
    };

    velocityIC()=velocitySolution();
    pressureIC()=pressureSolution();

    velocityBC().addAllBoundaryIDs(velocitySolution());
  }
};

// Navier-Stokes 2D
template<typename FluidStateImp>
class NavierStokes2DProblem:public BaseProblem<FluidStateImp,DirichletCondition>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef NavierStokes2DProblem<FluidStateType> ThisType;
  typedef BaseProblem<FluidStateType,DirichletCondition> BaseType;

  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::VelocityDomainType VelocityDomainType;
  typedef typename BaseType::VelocityRangeType VelocityRangeType;
  typedef typename BaseType::PressureDomainType PressureDomainType;
  typedef typename BaseType::PressureRangeType PressureRangeType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::pressureSolution;
  using BaseType::pressureIC;
  using BaseType::mu;

  NavierStokes2DProblem(FluidStateType& fluidState):
    BaseType(fluidState,true,true,"Navier-Stokes 2D")
  {
    velocitySolution().function()=[&](const VelocityDomainType& x,double t,const EntityType& entity)
    {
      const auto muValue(mu(entity));
      VelocityRangeType value(0.0);
      value[0]=-std::cos(M_PI*x[0])*std::sin(M_PI*x[1])*std::exp(-2.0*M_PI*M_PI*muValue*t);
      value[1]=std::sin(M_PI*x[0])*std::cos(M_PI*x[1])*std::exp(-2.0*M_PI*M_PI*muValue*t);
      return value;
    };

    pressureSolution().function()=[&](const PressureDomainType& x,double t,const EntityType& entity)
    {
      typename BaseType::PressureRangeType value(0.0);
      value[0]=-0.25*(std::cos(2.0*M_PI*x[0])+std::cos(2.0*M_PI*x[1]))*std::exp(-4.0*M_PI*M_PI*mu(entity)*t);
      return value;
    };

    velocityIC()=velocitySolution();
    pressureIC()=pressureSolution();

    velocityBC().addAllBoundaryIDs(velocitySolution());
  }
};

// rising bubble
template<typename FluidStateImp>
class RisingBubbleProblem:public BaseProblem<FluidStateImp,DirichletCondition,FreeSlipCondition>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef RisingBubbleProblem<FluidStateType> ThisType;
  typedef BaseProblem<FluidStateType,DirichletCondition,FreeSlipCondition> BaseType;

  typedef typename BaseType::VelocityFunctionType VelocityFunctionType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::VelocityDomainType VelocityDomainType;
  typedef typename BaseType::VelocityRangeType VelocityRangeType;

  using BaseType::velocityBC;
  using BaseType::velocityRHS;
  using BaseType::rho;
  using BaseType::worlddim;

  RisingBubbleProblem(FluidStateType& fluidState):
    BaseType(fluidState,true,false,"rising bubble")
  {
    velocityRHS().function()=[&](const VelocityDomainType& ,double ,const EntityType& entity)
    {
      VelocityRangeType value(0.0);
      value[worlddim-1]=-0.98;
      value*=rho(entity);
      return value;
    };

    // Dirichlet on top/bottom
    velocityBC().addBC(2,VelocityFunctionType());
    velocityBC().addBC(4,VelocityFunctionType());
    // free-slip on faces normal to x
    this->template getVelocityBC<1>().addBC(3);
    this->template getVelocityBC<1>().addBC(5);
    // free-slip on faces normal to y (only 3d)
    if(worlddim==3)
    {
      this->template getVelocityBC<1>().addBC(6);
      this->template getVelocityBC<1>().addBC(7);
    }
  }
};

// Navier-Stokes test 1
template<typename FluidStateImp>
class NavierStokesTest1Problem:public BaseProblem<FluidStateImp,DirichletCondition>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef NavierStokesTest1Problem<FluidStateType> ThisType;
  typedef BaseProblem<FluidStateType,DirichletCondition> BaseType;

  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::VelocityDomainType VelocityDomainType;
  typedef typename BaseType::VelocityDomainType VelocityRangeType;
  typedef typename BaseType::PressureFunctionType PressureFunctionType;
  typedef typename BaseType::PressureDomainType PressureDomainType;
  typedef typename BaseType::PressureRangeType PressureRangeType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::velocityRHS;
  using BaseType::pressureSolution;
  using BaseType::pressureIC;
  using BaseType::mu_;
  using BaseType::gamma;
  using BaseType::rho;
  using BaseType::worlddim;
  using BaseType::fluidstate_;

  NavierStokesTest1Problem(FluidStateType& fluidState):
    BaseType(fluidState,true,true,"Navier-Stokes test 1"),r0_(0.5),alpha_(0.15)
  {
    velocityRHS().function()=[&](const VelocityDomainType& x,double ,const EntityType& entity)
    {
      VelocityRangeType value(x);
      value*=(std::pow(alpha_,2)*rho(entity));
      return value;
    };

    velocitySolution().function()=[&](const VelocityDomainType& x,double ,const EntityType& )
    {
      auto value(x);
      value*=alpha_;
      return value;
    };

    pressureSolution().function()=[&](const PressureDomainType& x,double t,const EntityType& )
    {
      const auto rt(exactRadius(t));
      const auto coeff((static_cast<double>(worlddim-1)/rt)*gamma()+2.0*alpha_*mu_.delta());
      const auto indicatorValue(x.two_norm()<=rt?1.0:0.0);
      auto value(coeff*(indicatorValue-(std::pow(4.0/3.0,worlddim-2)*M_PI*std::pow(rt,worlddim))/(std::pow(2.0,worlddim))));
      return value;
    };

    velocityIC()=velocitySolution();
    pressureIC()=pressureSolution();

    velocityBC().addAllBoundaryIDs(velocitySolution());
  }

  double exactRadius(double t) const
  {
    return std::exp(alpha_*t)*r0_;
  }

  private:
  const double r0_;
  const double alpha_;
};

// Navier-Stokes test 2
template<typename FluidStateImp>
class NavierStokesTest2Problem:public BaseProblem<FluidStateImp,DirichletCondition>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef NavierStokesTest1Problem<FluidStateType> ThisType;
  typedef BaseProblem<FluidStateType,DirichletCondition> BaseType;

  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::VelocityDomainType VelocityDomainType;
  typedef typename BaseType::VelocityDomainType VelocityRangeType;
  typedef typename BaseType::PressureFunctionType PressureFunctionType;
  typedef typename BaseType::PressureDomainType PressureDomainType;
  typedef typename BaseType::PressureRangeType PressureRangeType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::velocityRHS;
  using BaseType::pressureSolution;
  using BaseType::pressureIC;
  using BaseType::mu_;
  using BaseType::gamma;
  using BaseType::rho;
  using BaseType::worlddim;
  using BaseType::fluidstate_;

  NavierStokesTest2Problem(FluidStateType& fluidState):
    BaseType(fluidState,true,true,"Navier-Stokes test 2"),r0_(0.5),alpha1_(0.15),alpha2_(0.15)
  {
    velocityRHS().function()=[&](const VelocityDomainType& x,double t,const EntityType& entity)
    {
      VelocityRangeType value(x);
      const auto rt2(std::pow(exactRadius(t),2));
      const auto x2(std::pow(x.two_norm(),2));
      const auto indicatorValue(fluidstate_.meshManager().bulkInnerIndicatorFunction().contains(entity)?0.0:1.0);
      value*=((std::pow(alpha1_,2)+indicatorValue*(-2.0*alpha1_*alpha2_*rt2+2.0*alpha1_*alpha2_*(2.0*x2-rt2)+
               std::pow(alpha2_,2)*(x2-rt2)*(3.0*x2-rt2)))*rho(entity)-4.0*(worlddim+1)*alpha2_*mu_.outerValue()*indicatorValue);
      return value;
    };

    velocitySolution().function()=[&](const VelocityDomainType& x,double t,const EntityType& )
    {
      auto value(x);
      const auto rt2(std::pow(exactRadius(t),2));
      const auto x2(std::pow(x.two_norm(),2));
      const auto indicatorValue(x2<rt2?0.0:1.0);
      value*=(alpha1_+alpha2_*indicatorValue*(x2-rt2));
      return value;
    };

    pressureSolution().function()=[&](const PressureDomainType& x,double t,const EntityType& )
    {
      const auto rt(exactRadius(t));
      const auto coeff((static_cast<double>(worlddim-1)/rt)*gamma()+4.0*alpha2_*mu_.outerValue()*std::pow(rt,2));
      const auto indicatorValue(x.two_norm()<=rt?1.0:0.0);
      auto value(coeff*(indicatorValue-(std::pow(4.0/3.0,worlddim-2)*M_PI*std::pow(rt,worlddim))/(std::pow(2.0,worlddim))));
      return value;
    };

    velocityIC()=velocitySolution();
    pressureIC()=pressureSolution();

    velocityBC().addAllBoundaryIDs(velocitySolution());
  }

  double exactRadius(double t) const
  {
    return std::exp(alpha1_*t)*r0_;
  }

  private:
  const double r0_;
  const double alpha1_;
  const double alpha2_;
};

// Navier-Stokes expanding bubble
template<typename FluidStateImp>
class NavierStokesExpandingBubbleProblem:public ExpandingBubbleProblem<FluidStateImp>
{
  public:
  typedef FluidStateImp FluidStateType;
  typedef NavierStokesExpandingBubbleProblem<FluidStateType> ThisType;
  typedef ExpandingBubbleProblem<FluidStateType> BaseType;

  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::VelocityDomainType VelocityDomainType;

  using BaseType::velocityRHS;
  using BaseType::alpha_;
  using BaseType::worlddim;

  NavierStokesExpandingBubbleProblem(FluidStateType& fluidState):
    BaseType(fluidState)
  {
    velocityRHS().function()=[&](const VelocityDomainType& x,double ,const EntityType& entity)
    {
      auto value(x);
      value*=-std::pow(alpha_,2)*(worlddim-1);
      value/=std::pow(x.two_norm(),2*worlddim);
      return value;
    };
  }
};

}
}

#endif // DUNE_FEM_PROBLEMS_HH
