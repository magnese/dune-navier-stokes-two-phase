#ifndef DUNE_FEM_PROBLEMS_HH
#define DUNE_FEM_PROBLEMS_HH

#include <cmath>
#include <tuple>
#include <ostream>
#include <iostream>
#include <string>
#include <functional>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>

#include "boundarycondition.hh"

namespace Dune
{
namespace Fem
{

// base problem
template<class VelocityDiscreteSpaceImp,class PressureDiscreteSpaceImp,class CoupledMeshManagerImp,
         template<class,class,class> class... VelocityBCImp>
class BaseProblem
{
  public:
  typedef VelocityDiscreteSpaceImp VelocityDiscreteSpaceType;
  typedef PressureDiscreteSpaceImp PressureDiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef BaseProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType,VelocityBCImp...> ThisType;

  typedef typename VelocityDiscreteSpaceType::EntityType EntityType;

  typedef typename VelocityDiscreteSpaceType::DomainType VelocityDomainType;
  typedef typename VelocityDiscreteSpaceType::RangeType VelocityRangeType;
  typedef std::function<VelocityRangeType(const VelocityDomainType&,const double&,const EntityType&)> VelocityFunctionType;
  typedef typename PressureDiscreteSpaceType::DomainType PressureDomainType;
  typedef typename PressureDiscreteSpaceType::RangeType PressureRangeType;
  typedef std::function<PressureRangeType(const PressureDomainType&,const double&,const EntityType&)> PressureFunctionType;
  typedef std::tuple<VelocityBCImp<VelocityDiscreteSpaceType,VelocityDiscreteSpaceType,CoupledMeshManagerType>...> VelocityBCsType;

  BaseProblem(CoupledMeshManagerType& meshManager,const bool& isTimeDependent,const bool& hasExactSolution,const std::string& name):
    meshmanager_(meshManager),
    velocitybcs_(VelocityBCImp<VelocityDiscreteSpaceType,VelocityDiscreteSpaceType,CoupledMeshManagerType>(meshmanager_)...),
    istimedependent_(isTimeDependent),
    hasexactsolution_(hasExactSolution),
    name_(name),
    muinner_(Parameter::getValidValue<double>("MuInner",1.0,[](double val){return val>0.0;})),
    muouter_(Parameter::getValidValue<double>("MuOuter",1.0,[](double val){return val>0.0;})),
    gamma_(Parameter::getValue<double>("Gamma",1.0)),
    rhoinner_(Parameter::getValidValue<double>("RhoInner",1.0,[](double val){return val>=0.0;})),
    rhoouter_(Parameter::getValidValue<double>("RhoOuter",1.0,[](double val){return val>=0.0;})),
    nulldensity_((rhoinner_==0.0&&rhoouter_==0.0)?true:false)
  {
    // init all the functions with the null function
    velocityRHS()=[&](const VelocityDomainType& ,const double& ,const EntityType& )
    {
      VelocityRangeType value(0.0);
      return value;
    };
    velocityIC()=[&](const VelocityDomainType& ,const double& ,const EntityType& )
    {
      VelocityRangeType value(0.0);
      return value;
    };
    velocitySolution()=[&](const VelocityDomainType& ,const double& ,const EntityType& )
    {
      VelocityRangeType value(0.0);
      return value;
    };
    pressureIC()=[&](const PressureDomainType& ,const double& ,const EntityType& )
    {
      PressureRangeType value(0.0);
      return value;
    };
    pressureSolution()=[&](const PressureDomainType& ,const double& ,const EntityType& )
    {
      PressureRangeType value(0.0);
      return value;
    };
  }

  inline const bool& isTimeDependent() const
  {
    return istimedependent_;
  }
  inline const bool& hasExactSolution() const
  {
    return hasexactsolution_;
  }

  inline VelocityFunctionType& velocityIC()
  {
    return std::get<1>(velocity_);
  }
  inline const VelocityFunctionType& velocityIC() const
  {
    return std::get<1>(velocity_);
  }
  template<class DF>
  inline void velocityIC(DF& df) const
  {
    if(istimedependent_)
      interpolateAnalyticalFunction(velocityIC(),df,0.0);
  }
  inline VelocityFunctionType& velocityRHS()
  {
    return std::get<2>(velocity_);
  }
  inline const VelocityFunctionType& velocityRHS() const
  {
    return std::get<2>(velocity_);
  }
  inline VelocityFunctionType& velocitySolution()
  {
    return std::get<0>(velocity_);
  }
  inline const VelocityFunctionType& velocitySolution() const
  {
    return std::get<0>(velocity_);
  }
  template<class DF>
  inline void velocitySolution(DF& df,const double& t) const
  {
    if(hasexactsolution_)
      interpolateAnalyticalFunction(velocitySolution(),df,t);
  }

  inline PressureFunctionType& pressureIC()
  {
    return std::get<1>(pressure_);
  }
  inline const PressureFunctionType& pressureIC() const
  {
    return std::get<1>(pressure_);
  }
  template<class DF>
  inline void pressureIC(DF& df) const
  {
    if(istimedependent_)
      interpolateAnalyticalFunction(pressureIC(),df,0.0);
  }
  inline PressureFunctionType& pressureSolution()
  {
    return std::get<0>(pressure_);
  }
  inline const PressureFunctionType& pressureSolution() const
  {
    return std::get<0>(pressure_);
  }
  template<class DF>
  inline void pressureSolution(DF& df,const double& t) const
  {
    if(hasexactsolution_)
      interpolateAnalyticalFunction(pressureSolution(),df,t);
  }

  inline double mu(const EntityType& entity) const
  {
    const auto& indicator(meshmanager_.bulkIndicatorFunction());
    return muinner_*indicator(entity)+muouter_*(1.0-indicator(entity));
  }
  inline double deltaMu() const
  {
    return muouter_-muinner_;
  }
  inline double gamma() const
  {
    return gamma_;
  }
  inline double rho(const EntityType& entity) const
  {
    const auto& indicator(meshmanager_.bulkIndicatorFunction());
    return rhoinner_*indicator(entity)+rhoouter_*(1.0-indicator(entity));
  }
  inline double deltaRho() const
  {
    return rhoouter_-rhoinner_;
  }
  inline bool isDensityNull() const
  {
    return nulldensity_;
  }

  void printInfo(std::ostream& s=std::cout) const
  {
    s<<"Problem : "<<name_<<std::endl;
    s<<"mu_inner = "<<muinner_<<std::endl;
    s<<"mu_outer = "<<muouter_<<std::endl;
    s<<"gamma = "<<gamma_<<std::endl;
    s<<"rho_inner = "<<rhoinner_<<std::endl;
    s<<"rho_outer = "<<rhoouter_<<std::endl;
    s<<std::endl;
  }

  protected:
  CoupledMeshManagerType& meshmanager_;
  VelocityBCsType velocitybcs_;
  const bool istimedependent_;
  const bool hasexactsolution_;
  const std::string name_;
  const double muinner_;
  const double muouter_;
  const double gamma_;
  const double rhoinner_;
  const double rhoouter_;
  std::tuple<VelocityFunctionType,VelocityFunctionType,VelocityFunctionType> velocity_;
  std::tuple<PressureFunctionType,PressureFunctionType> pressure_;
  static constexpr auto worlddim=CoupledMeshManagerType::BulkGridType::dimensionworld;
  static constexpr auto numbcs_=std::tuple_size<VelocityBCsType>::value;
  const bool nulldensity_;

  template<typename AF,typename DF>
  void interpolateAnalyticalFunction(const AF& f,DF& df,const double& t) const
  {
    // create local analytical function
    typedef typename DF::DiscreteFunctionSpaceType DiscreteSpaceType;
    typedef LocalAnalyticalFunctionBinder<DiscreteSpaceType> LocalAnalyticalFunctionType;
    LocalAnalyticalFunctionType localAnalyticalFunction(f);
    localAnalyticalFunction.initialize(t);
    // create local function adapter
    typedef LocalFunctionAdapter<LocalAnalyticalFunctionType> AdaptedFunctionType;
    AdaptedFunctionType fAdapted("adapted function",localAnalyticalFunction,df.gridPart(),1);
    // interpolate adpated function over df
    interpolate(fAdapted,df);
  }

  template<std::size_t N>
  auto getVelocityBC()->decltype(std::get<N>(velocitybcs_))
  {
    return std::get<N>(velocitybcs_);
  }

  struct ApplyBCToOperator
  {
    template<typename Tuple,typename... Args>
    ApplyBCToOperator(Tuple&& t,Args&&... args)
    {
      Caller<Tuple,numbcs_,Args...> caller(t,args...);
    }
    template<typename Tuple,std::size_t pos,typename... Args>
    struct Caller
    {
      inline Caller(Tuple&& t,Args&&... args)
      {
        std::get<pos-1>(t).applyToOperator(args...);
        Caller<Tuple,pos-1,Args...>(t,args...);
      }
    };
    template<typename Tuple,typename... Args>
    struct Caller<Tuple,0,Args...>
    {
      inline Caller(Tuple&& ,Args&&... )
      {}
    };
  };

  struct ApplyBCToRHS
  {
    template<typename Tuple,typename... Args>
    ApplyBCToRHS(Tuple&& t,Args&&... args)
    {
      Caller<Tuple,numbcs_,Args...> caller(t,args...);
    }
    template<typename Tuple,std::size_t pos,typename... Args>
    struct Caller
    {
      inline Caller(Tuple&& t,Args&&... args)
      {
        std::get<pos-1>(t).applyToRHS(args...);
        Caller<Tuple,pos-1,Args...>(t,args...);
      }
    };
    template<typename Tuple,typename... Args>
    struct Caller<Tuple,0,Args...>
    {
      inline Caller(Tuple&& ,Args&&... )
      {}
    };
  };

  public:
  auto velocityBC()->decltype(this->getVelocityBC<0>())
  {
    return this->getVelocityBC<0>();
  }

  template<typename... Args>
  inline void applyBCToOperator(Args&&... args)
  {
    ApplyBCToOperator apply(velocitybcs_,args...);
  }

  template<typename... Args>
  inline void applyBCToRHS(Args&&... args)
  {
    ApplyBCToRHS apply(velocitybcs_,args...);
  }
};

// Stokes test 1
template<class VelocityDiscreteSpaceImp,class PressureDiscreteSpaceImp,class CoupledMeshManagerImp>
class StokesTest1Problem:public BaseProblem<VelocityDiscreteSpaceImp,PressureDiscreteSpaceImp,CoupledMeshManagerImp,DirichletCondition>
{
  public:
  typedef VelocityDiscreteSpaceImp VelocityDiscreteSpaceType;
  typedef PressureDiscreteSpaceImp PressureDiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef StokesTest1Problem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType> ThisType;
  typedef BaseProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType,DirichletCondition> BaseType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::velocityRHS;
  using BaseType::mu;

  StokesTest1Problem(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,false,true,"Stokes test 1")
  {
    velocityRHS()=[&](const typename BaseType::VelocityDomainType& x,const double& ,const typename BaseType::EntityType& entity)
    {
      typename BaseType::VelocityRangeType value(0.0);
      value[1]=-6.0*mu(entity)*x[0];
      return value;
    };

    velocitySolution()=[&](const typename BaseType::VelocityDomainType& x,const double& ,const typename BaseType::EntityType& )
    {
      typename BaseType::VelocityRangeType value(0.0);
      value[1]=pow(x[0],3.0)-x[0];
      return value;
    };

    velocityIC()=velocitySolution();

    velocityBC().addBC(2,velocitySolution());
    velocityBC().addBC(3,velocitySolution());
    velocityBC().addBC(4,velocitySolution());
    velocityBC().addBC(5,velocitySolution());
  }
};

// Stokes test 2
template<class VelocityDiscreteSpaceImp,class PressureDiscreteSpaceImp,class CoupledMeshManagerImp>
class StokesTest2Problem:public BaseProblem<VelocityDiscreteSpaceImp,PressureDiscreteSpaceImp,CoupledMeshManagerImp,DirichletCondition>
{
  public:
  typedef VelocityDiscreteSpaceImp VelocityDiscreteSpaceType;
  typedef PressureDiscreteSpaceImp PressureDiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef StokesTest2Problem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType> ThisType;
  typedef BaseProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType,DirichletCondition> BaseType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::velocityRHS;
  using BaseType::mu;

  StokesTest2Problem(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,false,true,"Stokes test 2")
  {
    velocityRHS()=[&](const typename BaseType::VelocityDomainType& x,const double& ,const typename BaseType::EntityType& entity)
    {
      const auto muValue(mu(entity));
      typename BaseType::VelocityRangeType value(0.0);
      value[0]=-12.0*muValue*pow(x[1],2.0);
      value[1]=-6.0*muValue*x[0];
      return value;
    };

    velocitySolution()=[&](const typename BaseType::VelocityDomainType& x,const double& ,const typename BaseType::EntityType& )
    {
      typename BaseType::VelocityRangeType value(0.0);
      value[0]=1.0+pow(x[1],4.0);
      value[1]=pow(x[0],3.0)-x[0];
      return value;
    };

    velocityIC()=velocitySolution();

    velocityBC().addBC(2,velocitySolution());
    velocityBC().addBC(3,velocitySolution());
    velocityBC().addBC(4,velocitySolution());
    velocityBC().addBC(5,velocitySolution());
  }
};

// stationary bubble
template<class VelocityDiscreteSpaceImp,class PressureDiscreteSpaceImp,class CoupledMeshManagerImp>
class StationaryBubbleProblem:public BaseProblem<VelocityDiscreteSpaceImp,PressureDiscreteSpaceImp,CoupledMeshManagerImp,DirichletCondition>
{
  public:
  typedef VelocityDiscreteSpaceImp VelocityDiscreteSpaceType;
  typedef PressureDiscreteSpaceImp PressureDiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef StationaryBubbleProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType> ThisType;
  typedef BaseProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType,DirichletCondition> BaseType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::pressureSolution;
  using BaseType::pressureIC;
  using BaseType::gamma;
  using BaseType::worlddim;
  using BaseType::meshmanager_;

  StationaryBubbleProblem(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,true,true,"stationary bubble"),r0_(0.5)
  {
    pressureSolution()=[&](const typename BaseType::PressureDomainType& ,const double& ,const typename BaseType::EntityType& entity)
    {
      const auto rt(exactRadius());
      const auto lambda(gamma()*static_cast<double>(worlddim-1)/rt);
      const auto& indicator(meshmanager_.bulkIndicatorFunction());
      auto value(lambda*(indicator(entity)-pow(4.0/3.0,worlddim-2)*M_PI*pow(rt,worlddim)*pow(2,-worlddim)));
      return value;
    };

    pressureIC()=pressureSolution();

    velocityBC().addBC(2,velocitySolution());
    velocityBC().addBC(3,velocitySolution());
    velocityBC().addBC(4,velocitySolution());
    velocityBC().addBC(5,velocitySolution());
    velocityBC().addBC(6,velocitySolution());
    velocityBC().addBC(7,velocitySolution());
  }

  template<typename... Args>
  inline double exactRadius(const Args&... ) const
  {
    return 0.5;
  }

  private:
  const double r0_;
};

// expanding bubble
template<class VelocityDiscreteSpaceImp,class PressureDiscreteSpaceImp,class CoupledMeshManagerImp>
class ExpandingBubbleProblem:public BaseProblem<VelocityDiscreteSpaceImp,PressureDiscreteSpaceImp,CoupledMeshManagerImp,DirichletCondition>
{
  public:
  typedef VelocityDiscreteSpaceImp VelocityDiscreteSpaceType;
  typedef PressureDiscreteSpaceImp PressureDiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef ExpandingBubbleProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType> ThisType;
  typedef BaseProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType,DirichletCondition> BaseType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::pressureSolution;
  using BaseType::pressureIC;
  using BaseType::deltaMu;
  using BaseType::gamma;
  using BaseType::worlddim;
  using BaseType::meshmanager_;

  ExpandingBubbleProblem(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,true,true,"expanding bubble"),r0_(0.5),alpha_(0.15)
  {
    velocitySolution()=[&](const typename BaseType::VelocityDomainType& x,const double& ,const typename BaseType::EntityType& )
    {
      auto value(x);
      value*=alpha_;
      value/=pow(x.two_norm(),worlddim);
      return value;
    };

    pressureSolution()=[&](const typename BaseType::PressureDomainType& ,const double& t,const typename BaseType::EntityType& entity)
    {
      const auto rt(exactRadius(t));
      const auto lambda(static_cast<double>(worlddim-1)*(gamma()/rt+2*alpha_*deltaMu()*pow(rt,-worlddim)));
      const auto& indicator(meshmanager_.bulkIndicatorFunction());
      auto value(lambda*(indicator(entity)-
                         (pow(4.0/3.0,worlddim-2)*M_PI*pow(rt,worlddim)-pow(2.0/3.0,worlddim))/(pow(2,worlddim)-pow(2.0/3.0,worlddim))));
      return value;
    };

    velocityIC()=velocitySolution();
    pressureIC()=pressureSolution();

    velocityBC().addBC(2,velocitySolution());
    velocityBC().addBC(3,velocitySolution());
    velocityBC().addBC(4,velocitySolution());
    velocityBC().addBC(5,velocitySolution());
    velocityBC().addBC(6,velocitySolution());
    velocityBC().addBC(7,velocitySolution());
    velocityBC().addBC(8,velocitySolution());
    velocityBC().addBC(9,velocitySolution());
  }

  inline double exactRadius(const double& t) const
  {
    return pow(pow(r0_,worlddim)+alpha_*t*static_cast<double>(worlddim),1.0/static_cast<double>(worlddim));
  }

  private:
  const double r0_;
  const double alpha_;
};

// shear flow
template<class VelocityDiscreteSpaceImp,class PressureDiscreteSpaceImp,class CoupledMeshManagerImp>
class ShearFlowProblem:public BaseProblem<VelocityDiscreteSpaceImp,PressureDiscreteSpaceImp,CoupledMeshManagerImp,DirichletCondition>
{
  public:
  typedef VelocityDiscreteSpaceImp VelocityDiscreteSpaceType;
  typedef PressureDiscreteSpaceImp PressureDiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef ShearFlowProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType> ThisType;
  typedef BaseProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType,DirichletCondition> BaseType;

  using BaseType::velocityBC;
  using BaseType::worlddim;

  ShearFlowProblem(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,true,false,"shear flow")
  {
    f_=[&](const typename BaseType::VelocityDomainType& x,const double& ,const typename BaseType::EntityType& )
    {
      typename BaseType::VelocityDomainType value(0.0);
      value[0]=x[worlddim-1];
      return value;
    };

    velocityBC().addBC(2,f_);
    velocityBC().addBC(3,f_);
    velocityBC().addBC(4,f_);
    velocityBC().addBC(5,f_);
    velocityBC().addBC(6,f_);
    velocityBC().addBC(7,f_);
  }

  private:
  typename BaseType::VelocityFunctionType f_;
};

// stationary Navier-Stokes
template<class VelocityDiscreteSpaceImp,class PressureDiscreteSpaceImp,class CoupledMeshManagerImp>
class StationaryNavierStokesProblem:public BaseProblem<VelocityDiscreteSpaceImp,PressureDiscreteSpaceImp,CoupledMeshManagerImp,
                                                       DirichletCondition>
{
  public:
  typedef VelocityDiscreteSpaceImp VelocityDiscreteSpaceType;
  typedef PressureDiscreteSpaceImp PressureDiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef StationaryNavierStokesProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType> ThisType;
  typedef BaseProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType,DirichletCondition> BaseType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::velocityRHS;
  using BaseType::pressureSolution;
  using BaseType::pressureIC;
  using BaseType::mu;

  StationaryNavierStokesProblem(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,false,true,"stationary Navier-Stokes")
  {
    velocityRHS()=[&](const typename BaseType::VelocityDomainType& x,const double& ,const typename BaseType::EntityType& entity)
    {
      const auto muValue(mu(entity));
      typename BaseType::VelocityRangeType value(0.0);
      value[0]=-2.0*pow(M_PI,2.0)*muValue*cos(M_PI*x[0])*sin(M_PI*x[1])+M_PI/2.0*sin(2.0*M_PI*x[0]);
      value[1]=2.0*pow(M_PI,2.0)*muValue*sin(M_PI*x[0])*cos(M_PI*x[1])+M_PI/2.0*sin(2.0*M_PI*x[1]);
      return value;
    };

    velocitySolution()=[&](const typename BaseType::VelocityDomainType& x,const double& ,const typename BaseType::EntityType& )
    {
      typename BaseType::VelocityRangeType value(0.0);
      value[0]=-cos(M_PI*x[0])*sin(M_PI*x[1]);
      value[1]=sin(M_PI*x[0])*cos(M_PI*x[1]);
      return value;
    };

    pressureSolution()=[&](const typename BaseType::PressureDomainType& x,const double& ,const typename BaseType::EntityType& )
    {
      typename BaseType::PressureRangeType value(0.0);
      value[0]=-0.25*(cos(2.0*M_PI*x[0])+cos(2.0*M_PI*x[1]));
      return value;
    };

    velocityIC()=velocitySolution();
    pressureIC()=pressureSolution();

    velocityBC().addBC(2,velocitySolution());
    velocityBC().addBC(3,velocitySolution());
    velocityBC().addBC(4,velocitySolution());
    velocityBC().addBC(5,velocitySolution());
  }
};

// Navier-Stokes 2D
template<class VelocityDiscreteSpaceImp,class PressureDiscreteSpaceImp,class CoupledMeshManagerImp>
class NavierStokes2DProblem:public BaseProblem<VelocityDiscreteSpaceImp,PressureDiscreteSpaceImp,CoupledMeshManagerImp,DirichletCondition>
{
  public:
  typedef VelocityDiscreteSpaceImp VelocityDiscreteSpaceType;
  typedef PressureDiscreteSpaceImp PressureDiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef NavierStokes2DProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType> ThisType;
  typedef BaseProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType,DirichletCondition> BaseType;

  using BaseType::velocitySolution;
  using BaseType::velocityBC;
  using BaseType::velocityIC;
  using BaseType::pressureSolution;
  using BaseType::pressureIC;
  using BaseType::mu;

  NavierStokes2DProblem(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,true,true,"Navier-Stokes 2D")
  {
    velocitySolution()=[&](const typename BaseType::VelocityDomainType& x,const double& t,const typename BaseType::EntityType& entity)
    {
      const auto muValue(mu(entity));
      typename BaseType::VelocityRangeType value(0.0);
      value[0]=-cos(M_PI*x[0])*sin(M_PI*x[1])*exp(-2.0*M_PI*M_PI*muValue*t);
      value[1]=sin(M_PI*x[0])*cos(M_PI*x[1])*exp(-2.0*M_PI*M_PI*muValue*t);
      return value;
    };

    pressureSolution()=[&](const typename BaseType::PressureDomainType& x,const double& t,const typename BaseType::EntityType& entity)
    {
      typename BaseType::PressureRangeType value(0.0);
      value[0]=-0.25*(cos(2.0*M_PI*x[0])+cos(2.0*M_PI*x[1]))*exp(-4.0*M_PI*M_PI*mu(entity)*t);
      return value;
    };

    velocityIC()=velocitySolution();
    pressureIC()=pressureSolution();

    velocityBC().addBC(2,velocitySolution());
    velocityBC().addBC(3,velocitySolution());
    velocityBC().addBC(4,velocitySolution());
    velocityBC().addBC(5,velocitySolution());
  }
};

// rising bubble 2D
template<class VelocityDiscreteSpaceImp,class PressureDiscreteSpaceImp,class CoupledMeshManagerImp>
class RisingBubble2DProblem:public BaseProblem<VelocityDiscreteSpaceImp,PressureDiscreteSpaceImp,CoupledMeshManagerImp,DirichletCondition,
                                               FreeSlipCondition>
{
  public:
  typedef VelocityDiscreteSpaceImp VelocityDiscreteSpaceType;
  typedef PressureDiscreteSpaceImp PressureDiscreteSpaceType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef RisingBubble2DProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType> ThisType;
  typedef BaseProblem<VelocityDiscreteSpaceType,PressureDiscreteSpaceType,CoupledMeshManagerType,DirichletCondition,FreeSlipCondition>
    BaseType;

  using BaseType::velocityBC;
  using BaseType::velocityRHS;
  using BaseType::rho;

  RisingBubble2DProblem(CoupledMeshManagerType& meshManager):
    BaseType(meshManager,true,false,"rising bubble 2D")
  {
    typedef typename BaseType::VelocityDomainType VelocityDomainType;
    typedef typename BaseType::VelocityRangeType VelocityRangeType;
    typedef typename BaseType::EntityType EntityType;

    velocityRHS()=[&](const VelocityDomainType& ,const double& ,const EntityType& entity)
    {
      typename BaseType::VelocityRangeType value(-0.98);
      value*=rho(entity);
      return value;
    };

    velocityBC().addBC(2,[&](const VelocityDomainType& ,const double& ,const EntityType& )
                            {
                              return VelocityRangeType(0.0);
                            });
    velocityBC().addBC(4,[&](const VelocityDomainType& ,const double& ,const EntityType& )
                            {
                              return VelocityRangeType(0.0);
                            });
    this->template getVelocityBC<1>().addBC(3,[&](const VelocityDomainType& ,const double& ,const EntityType& )
                                                 {
                                                   VelocityRangeType value(1.0);
                                                   value[0]=0.0;
                                                   return value;
                                                 });
    this->template getVelocityBC<1>().addBC(5,[&](const VelocityDomainType& ,const double& ,const EntityType& )
                                                 {
                                                   VelocityRangeType value(1.0);
                                                   value[0]=0.0;
                                                   return value;
                                                 });
  }
};

}
}

#endif // DUNE_FEM_PROBLEMS_HH
