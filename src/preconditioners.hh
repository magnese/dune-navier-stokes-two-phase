#ifndef DUNE_FEM_PRECONDITIONERS_HH
#define DUNE_FEM_PRECONDITIONERS_HH

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/fem/solver/umfpacksolver.hh>

#include <vector>

namespace Dune
{
namespace Fem
{

template<class Oper11T,class Oper12T,class Oper22T>
class StokesPrecond:public Dune::Preconditioner<
  typename ISTLBlockVectorDiscreteFunction<TupleDiscreteFunctionSpace<typename Oper11T::DomainSpaceType,
                                                                      typename Oper12T::DomainSpaceType>>::DofStorageType,
  typename ISTLBlockVectorDiscreteFunction<TupleDiscreteFunctionSpace<typename Oper11T::DomainSpaceType,
                                                                      typename Oper12T::DomainSpaceType>>::DofStorageType>
{
  public:
  typedef Oper11T Oper11Type;
  typedef Oper12T Oper12Type;
  typedef Oper22T Oper22Type;

  typedef typename Oper11Type::DomainFunctionType DomainFunction1Type;
  typedef typename DomainFunction1Type::DiscreteFunctionSpaceType DomainSpace1Type;
  typedef typename Oper22Type::DomainFunctionType DomainFunction2Type;
  typedef typename DomainFunction2Type::DiscreteFunctionSpaceType DomainSpace2Type;
  typedef TupleDiscreteFunctionSpace<DomainSpace1Type,DomainSpace2Type> CombinedDiscreteFunctionSpaceType;
  typedef ISTLBlockVectorDiscreteFunction<CombinedDiscreteFunctionSpaceType> DomainFunctionType;
  typedef DomainFunctionType RangeFunctionType;

  typedef typename DomainFunctionType::DofStorageType domain_type;
  typedef typename RangeFunctionType::DofStorageType range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};

  // constructor
  explicit inline StokesPrecond(const Oper11Type& op11,const Oper12Type& op12,const Oper22Type& op22):
    op11_(op11),op12_(op12),op22_(op22),invop11_(op11),invop22_(op22),d1_(op11_.domainSpace().size(),0),x1_(op11_.domainSpace().size(),0),
    d2_(op12_.domainSpace().size(),0),x2_(op12_.domainSpace().size(),0)
  {}

  inline ~StokesPrecond()
  {}

  inline void pre(domain_type& ,range_type& )
  {
    invop11_.prepare();
    invop22_.prepare();
  }

  void apply(domain_type& x,const range_type& d)
  {
    // extract pointers to the vectors
    const double* d1Ptr(d1_.data());
    double* x1Ptr(x1_.data());
    const double* d2Ptr(d2_.data());
    double* x2Ptr(x2_.data());

    // split d into d1 and d2
    std::size_t pos(0);
    for(auto it=d1_.begin();it!=d1_.end();++it,++pos)
      (*it)=d[pos];
    for(auto it=d2_.begin();it!=d2_.end();++it,++pos)
      (*it)=d[pos];

    // solve -op22*x2 = d2 and copy x2 into x
    invop22_.apply(d2Ptr,x2Ptr);
    pos=op11_.domainSpace().size();
    for(auto it=x2_.begin();it!=x2_.end();++it,++pos)
      x[pos]=(-1.0)*(*it);

    // d1 = d1 - op12*x2
    op12_.systemMatrix().multOEM(x2Ptr,x1Ptr);
    auto x1It(x1_.begin());
    for(auto it=d1_.begin();it!=d1_.end();++it,++x1It)
      (*it)+=(*x1It);

    // solve op11*x1 = d1 and copy x1 into x
    invop11_.apply(d1Ptr,x1Ptr);
    pos=0;
    for(auto it=x1_.begin();it!=x1_.end();++it,++pos)
      x[pos]=(*it);
  }

  inline void post(domain_type& )
  {
    invop11_.finalize();
    invop22_.finalize();
  }

  private:
  const Oper11Type& op11_;
  const Oper12Type& op12_;
  const Oper22Type& op22_;
  UMFPACKOp<DomainFunction1Type,Oper11Type> invop11_;
  UMFPACKOp<DomainFunction2Type,Oper22Type> invop22_;
  std::vector<double> d1_;
  std::vector<double> x1_;
  std::vector<double> d2_;
  std::vector<double> x2_;
};


template<class Oper11T,class Oper12T,class Oper22T,class Oper13T,class Oper33T>
class ExtendedStokesPrecond:public Dune::Preconditioner<
  typename ISTLBlockVectorDiscreteFunction<TupleDiscreteFunctionSpace<typename Oper11T::DomainSpaceType,
                                                                      typename Oper12T::DomainSpaceType,
                                                                      typename Oper13T::DomainSpaceType>>::DofStorageType,
  typename ISTLBlockVectorDiscreteFunction<TupleDiscreteFunctionSpace<typename Oper11T::DomainSpaceType,
                                                                      typename Oper12T::DomainSpaceType,
                                                                      typename Oper13T::DomainSpaceType>>::DofStorageType>
{
  public:
  typedef Oper11T Oper11Type;
  typedef Oper12T Oper12Type;
  typedef Oper13T Oper13Type;
  typedef Oper22T Oper22Type;
  typedef Oper33T Oper33Type;

  typedef typename Oper11Type::DomainFunctionType DomainFunction1Type;
  typedef typename DomainFunction1Type::DiscreteFunctionSpaceType DomainSpace1Type;
  typedef typename Oper12Type::DomainFunctionType DomainFunction2Type;
  typedef typename DomainFunction2Type::DiscreteFunctionSpaceType DomainSpace2Type;
  typedef typename Oper13Type::DomainFunctionType DomainFunction3Type;
  typedef typename DomainFunction3Type::DiscreteFunctionSpaceType DomainSpace3Type;
  typedef TupleDiscreteFunctionSpace<DomainSpace1Type,DomainSpace2Type,DomainSpace3Type> CombinedDiscreteFunctionSpaceType;
  typedef ISTLBlockVectorDiscreteFunction<CombinedDiscreteFunctionSpaceType> DomainFunctionType;
  typedef DomainFunctionType RangeFunctionType;

  typedef typename DomainFunctionType::DofStorageType domain_type;
  typedef typename RangeFunctionType::DofStorageType range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};

  // constructor
  explicit inline ExtendedStokesPrecond(const Oper11Type& op11,const Oper12Type& op12,const Oper22Type& op22,const Oper13Type& op13,
                                        const Oper33Type& op33):
    op11_(op11),op12_(op12),op22_(op22),op13_(op13),op33_(op33),invop11_(op11),invop22_(op22),invop33_(op33),
    d1_(op11_.domainSpace().size(),0),x1_(op11_.domainSpace().size(),0),d2_(op12_.domainSpace().size(),0),x2_(op12_.domainSpace().size(),0),
    d3_(op13_.domainSpace().size(),0),x3_(op13_.domainSpace().size(),0)
  {}

  inline ~ExtendedStokesPrecond()
  {}

  inline void pre(domain_type& ,range_type& )
  {
    invop11_.prepare();
    invop22_.prepare();
    invop33_.prepare();
  }

  void apply(domain_type& x,const range_type& d)
  {
    // extract pointers to the vectors
    const double* d1Ptr(d1_.data());
    double* x1Ptr(x1_.data());
    const double* d2Ptr(d2_.data());
    double* x2Ptr(x2_.data());
    const double* d3Ptr(d3_.data());
    double* x3Ptr(x3_.data());

    // split d into d1, d2 and d3
    std::size_t pos(0);
    for(auto it=d1_.begin();it!=d1_.end();++it,++pos)
      (*it)=d[pos];
    for(auto it=d2_.begin();it!=d2_.end();++it,++pos)
      (*it)=d[pos];
    for(auto it=d3_.begin();it!=d3_.end();++it,++pos)
      (*it)=d[pos];

    // solve -op22*x2 = d2 and copy x2 into x
    invop22_.apply(d2Ptr,x2Ptr);
    pos=op11_.domainSpace().size();
    for(auto it=x2_.begin();it!=x2_.end();++it,++pos)
      x[pos]=(-1.0)*(*it);

    // solve -op33*x3 = d3 and copy x3 into x
    invop33_.apply(d3Ptr,x3Ptr);
    for(auto it=x3_.begin();it!=x3_.end();++it,++pos)
      x[pos]=(-1.0)*(*it);

    // d1 = d1 - op12*x2 - op13*x3
    op12_.systemMatrix().multOEM(x2Ptr,x1Ptr);
    auto x1It(x1_.begin());
    for(auto it=d1_.begin();it!=d1_.end();++it,++x1It)
      (*it)+=(*x1It);
    op13_.systemMatrix().multOEM(x3Ptr,x1Ptr);
    x1It=x1_.begin();
    for(auto it=d1_.begin();it!=d1_.end();++it,++x1It)
      (*it)+=(*x1It);

    // solve op11*x1 = d1 and copy x1 into x
    invop11_.apply(d1Ptr,x1Ptr);
    pos=0;
    for(auto it=x1_.begin();it!=x1_.end();++it,++pos)
      x[pos]=(*it);
  }

  inline void post(domain_type& )
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
  UMFPACKOp<DomainFunction1Type,Oper11Type> invop11_;
  UMFPACKOp<DomainFunction2Type,Oper22Type> invop22_;
  UMFPACKOp<DomainFunction3Type,Oper33Type> invop33_;
  std::vector<double> d1_;
  std::vector<double> x1_;
  std::vector<double> d2_;
  std::vector<double> x2_;
  std::vector<double> d3_;
  std::vector<double> x3_;
};


template<class Oper11T,class Oper12T,class Oper21T,class Oper22T>
class IdPrecond:public Dune::Preconditioner<
  typename ISTLBlockVectorDiscreteFunction<TupleDiscreteFunctionSpace<typename Oper11T::DomainSpaceType,
                                                                      typename Oper12T::DomainSpaceType>>::DofStorageType,
  typename ISTLBlockVectorDiscreteFunction<TupleDiscreteFunctionSpace<typename Oper11T::DomainSpaceType,
                                                                      typename Oper12T::DomainSpaceType>>::DofStorageType>
{
  public:
  typedef Oper11T Oper11Type;
  typedef Oper12T Oper12Type;
  typedef Oper21T Oper21Type;
  typedef Oper22T Oper22Type;

  typedef typename Oper11Type::DomainFunctionType DomainFunction1Type;
  typedef typename DomainFunction1Type::DiscreteFunctionSpaceType DomainSpace1Type;
  typedef typename Oper22Type::DomainFunctionType DomainFunction2Type;
  typedef typename DomainFunction2Type::DiscreteFunctionSpaceType DomainSpace2Type;
  typedef TupleDiscreteFunctionSpace<DomainSpace1Type,DomainSpace2Type> CombinedDiscreteFunctionSpaceType;
  typedef ISTLBlockVectorDiscreteFunction<CombinedDiscreteFunctionSpaceType> DomainFunctionType;
  typedef DomainFunctionType RangeFunctionType;

  typedef typename DomainFunctionType::DofStorageType domain_type;
  typedef typename RangeFunctionType::DofStorageType range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};

  // constructor
  explicit inline IdPrecond(const Oper11Type& ,const Oper12Type& ,const Oper21Type& ,const Oper22Type& )
  {}

  inline void pre(domain_type& ,range_type& )
  {}

  inline void apply(domain_type& x,const range_type& d)
  {
    for(auto i=0;i!=d.size();++i)
      x[i]=d[i];
  }

  inline void post(domain_type& )
  {}
};

}
}

#endif //DUNE_FEM_PRECONDITIONERS_HH
