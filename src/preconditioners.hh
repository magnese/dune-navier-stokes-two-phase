#ifndef DUNE_FEM_PRECONDITIONERS_HH
#define DUNE_FEM_PRECONDITIONERS_HH

#include <dune/common/hybridutilities.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/fem/solver/spqrsolver.hh>
#include <dune/fem/solver/umfpacksolver.hh>

#include <tuple>
#include <type_traits>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DFT,typename Oper11T,typename Oper12T,typename Oper22T>
class StokesPrecond:public Dune::Preconditioner<DFT,DFT>
{
  public:
  typedef Oper11T Oper11Type;
  typedef Oper12T Oper12Type;
  typedef Oper22T Oper22Type;

  typedef DFT domain_type;
  typedef domain_type range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};

  StokesPrecond(const Oper11Type& op11,const Oper12Type& op12,const Oper22Type& op22):
    op11_(op11),op12_(op12),op22_(op22),invop11_(op11),invop22_(op22)
  {}

  void pre(domain_type& ,range_type& )
  {
    invop11_.prepare();
    invop22_.prepare();
  }

  void apply(domain_type& x,const range_type& b)
  {
    // get reference sub-df's
    auto& x1(x.template subDiscreteFunction<0>());
    auto& x2(x.template subDiscreteFunction<1>());
    const auto& b1(b.template subDiscreteFunction<0>());
    const auto& b2(b.template subDiscreteFunction<1>());

    // x2 = - op22^-1 * b2
    invop22_.apply(b2,x2);
    x2*=-1.0;

    // temp1 = b1 - op12 * x2
    op12_(x2,x1);
    auto temp1(b1);
    temp1-=x1;

    // x1 = op11^-1 * temp1
    invop11_.apply(temp1,x1);
  }

  void post(domain_type& )
  {
    invop11_.finalize();
    invop22_.finalize();
  }

  private:
  const Oper11Type& op11_;
  const Oper12Type& op12_;
  const Oper22Type& op22_;
  UMFPACKOp<typename Oper11Type::DomainFunctionType,Oper11Type> invop11_;
  UMFPACKOp<typename Oper22Type::DomainFunctionType,Oper22Type> invop22_;
};


template<typename DFT,typename Oper11T,typename Oper12T,typename Oper22T,typename Oper13T,typename Oper33T>
class ExtendedStokesPrecond:public Dune::Preconditioner<DFT,DFT>
{
  public:
  typedef Oper11T Oper11Type;
  typedef Oper12T Oper12Type;
  typedef Oper13T Oper13Type;
  typedef Oper22T Oper22Type;
  typedef Oper33T Oper33Type;

  typedef DFT domain_type;
  typedef domain_type range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};

  explicit ExtendedStokesPrecond(const Oper11Type& op11,const Oper12Type& op12,const Oper22Type& op22,const Oper13Type& op13,
                                 const Oper33Type& op33):
    op11_(op11),op12_(op12),op22_(op22),op13_(op13),op33_(op33),invop11_(op11),invop22_(op22),invop33_(op33)
  {}

  void pre(domain_type& ,range_type& )
  {
    invop11_.prepare();
    invop22_.prepare();
    invop33_.prepare();
  }

  void apply(domain_type& x,const range_type& b)
  {
    // get reference sub-df's
    auto& x1(x.template subDiscreteFunction<0>());
    auto& x2(x.template subDiscreteFunction<1>());
    auto& x3(x.template subDiscreteFunction<2>());
    const auto& b1(b.template subDiscreteFunction<0>());
    const auto& b2(b.template subDiscreteFunction<1>());
    const auto& b3(b.template subDiscreteFunction<2>());

    // x2 = - op22^-1 * b2
    invop22_.apply(b2,x2);
    x2*=-1.0;

    // x3 = - op33^-1 * b3
    invop33_.apply(b3,x3);
    x3*=-1.0;

    // b1 = b1 - op12*x2 - op13*x3
    op12_(x2,x1);
    auto temp1(b1);
    temp1-=x1;
    op13_(x3,x1);
    temp1-=x1;

    // x1 = op11^-1 * temp1
    invop11_.apply(temp1,x1);
  }

  void post(domain_type& )
  {
    invop11_.finalize();
    invop22_.finalize();
    invop33_.finalize();
  }

  private:
  const Oper11Type& op11_;
  const Oper12Type& op12_;
  const Oper22Type& op22_;
  const Oper13Type& op13_;
  const Oper33Type& op33_;
  UMFPACKOp<typename Oper11Type::DomainFunctionType,Oper11Type> invop11_;
  UMFPACKOp<typename Oper22Type::DomainFunctionType,Oper22Type> invop22_;
  UMFPACKOp<typename Oper33Type::DomainFunctionType,Oper33Type> invop33_;
};


template<typename OperT,template<typename ,typename ,bool > typename InvOperT=UMFPACKOp>
class DirectPrecond:public Dune::Preconditioner<typename OperT::DiscreteFunctionType,typename OperT::DiscreteFunctionType>
{
  public:
  typedef OperT OperType;
  typedef typename OperType::DiscreteFunctionType domain_type;
  typedef domain_type range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};
  typedef InvOperT<domain_type,OperType,false> InvOperType;

  explicit DirectPrecond(OperType& op):
    op_(op),invop_(op_),usedoctoring_((std::is_same<InvOperType,SPQROp<domain_type,OperType,false>>::value)?false:true),
    b_(op_.domainSpace().size(),0),x_(op_.domainSpace().size(),0)
  {}

  void pre(domain_type& ,range_type& )
  {
    if(usedoctoring_)
      op_.applyDoctoring();
    invop_.prepare();
  }

  void apply(domain_type& x,const range_type& b)
  {
    // copy b into b_
    auto bIt(b_.begin());
    Hybrid::forEach(typename range_type::Sequence{},[&](auto i){bIt=std::copy(std::get<i>(b).dbegin(),std::get<i>(b).dend(),bIt);});

    // solve system
    invop_.apply(b_.data(),x_.data());

    // copy x_ into x
    auto xIt(x_.begin());
    Hybrid::forEach(typename domain_type::Sequence{},[&](auto i){for(auto& dof:dofs(std::get<i>(x))) dof=(*(xIt++));});
  }

  void post(domain_type& )
  {
    invop_.finalize();
  }

  private:
  OperType& op_;
  InvOperType invop_;
  const bool usedoctoring_;
  std::vector<double> b_;
  std::vector<double> x_;
};


template<typename DFT>
class IdPrecond:public Dune::Preconditioner<DFT,DFT>
{
  public:
  typedef DFT domain_type;
  typedef domain_type range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};

  IdPrecond() = default;

  void pre(domain_type& ,range_type& )
  {}

  void apply(domain_type& x,const range_type& b)
  {
    x=b;
  }

  void post(domain_type& )
  {}
};

}
}

#endif // DUNE_FEM_PRECONDITIONERS_HH
