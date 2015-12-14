#ifndef DUNE_FEM_OPERATORWRAPPER_HH
#define DUNE_FEM_OPERATORWRAPPER_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/istl/operators.hh>

#include "nulloperator.hh"

namespace Dune
{
namespace Fem
{

template<typename Oper11T,typename Oper12T,typename Oper21T,typename Oper22T>
class OperatorWrapper:public Dune::LinearOperator<
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
  typedef typename DomainFunctionType::RangeFieldType DomainFieldType;
  typedef typename RangeFunctionType::RangeFieldType RangeFieldType;

  typedef typename DomainFunctionType::DofStorageType domain_type;
  typedef typename RangeFunctionType::DofStorageType range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};

  typedef typename Oper22Type::RangeFunctionType RangeFunction2Type;
  typedef typename Oper11Type::RangeFunctionType RangeFunction1Type;

  // constructor
  OperatorWrapper(const Oper11Type& oper11,const Oper12Type& oper12,const Oper21Type& oper21,const Oper22Type& oper22):
    oper11_(oper11),oper12_(oper12),oper21_(oper21),oper22_(oper22)
  {}

  void apply(const domain_type& x,range_type& b) const
  {
    // split x into x1 and x2
    DomainFunction1Type x1("x1",oper11_.domainSpace());
    DomainFunction2Type x2("x2",oper22_.domainSpace());
    std::size_t counter(0);
    for(auto& dof:dofs(x1))
    {
      dof=x[counter];
      ++counter;
    }
    for(auto& dof:dofs(x2))
    {
      dof=x[counter];
      ++counter;
    }

    // create b1, b1Temp, b2 and b2Temp
    RangeFunction1Type b1("b1",oper11_.rangeSpace());
    RangeFunction1Type b1Temp("b1Temp",oper11_.rangeSpace());
    RangeFunction2Type b2("b2",oper22_.rangeSpace());
    RangeFunction2Type b2Temp("b2Temp",oper22_.rangeSpace());

    // b1 = A11*x1 + A12*x2
    if(isNullOperator<Oper11Type>::value)
      b1.clear();
    else
      oper11_(x1,b1);
    if(isNullOperator<Oper12Type>::value)
      b1Temp.clear();
    else
      oper12_(x2,b1Temp);
    b1+=b1Temp;

    // b2 = A21*x1 + A22*x2
    if(isNullOperator<Oper21Type>::value)
      b2.clear();
    else
      oper21_(x1,b2);
    if(isNullOperator<Oper22Type>::value)
      b2Temp.clear();
    else
      oper22_(x2,b2Temp);
    b2+=b2Temp;

    // merge b1 and b2 into b
    counter=0;
    for(auto& dof:dofs(b1))
    {
      b[counter]=dof;
      ++counter;
    }
    for(auto& dof:dofs(b2))
    {
      b[counter]=dof;
      ++counter;
    }
  }

  void applyscaleadd(field_type alpha,const domain_type& x,range_type& b) const
  {
    range_type temp(b);
    apply(x,temp);
    const auto size(oper11_.domainSpace().size()+oper22_.domainSpace().size());
    for(auto i=0;i!=size;++i)
      b[i]+=alpha*temp[i];
  }

  private:
  const Oper11Type& oper11_;
  const Oper12Type& oper12_;
  const Oper21Type& oper21_;
  const Oper22Type& oper22_;
};

template<typename Oper11T,typename Oper12T,typename Oper21T,typename Oper13T,typename Oper31T>
class ExtendedOperatorWrapper:public Dune::LinearOperator<
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
  typedef Oper21T Oper21Type;
  typedef Oper13T Oper13Type;
  typedef Oper31T Oper31Type;

  typedef typename Oper11Type::DomainFunctionType DomainFunction1Type;
  typedef typename DomainFunction1Type::DiscreteFunctionSpaceType DomainSpace1Type;
  typedef typename Oper12Type::DomainFunctionType DomainFunction2Type;
  typedef typename DomainFunction2Type::DiscreteFunctionSpaceType DomainSpace2Type;
  typedef typename Oper13Type::DomainFunctionType DomainFunction3Type;
  typedef typename DomainFunction3Type::DiscreteFunctionSpaceType DomainSpace3Type;
  typedef TupleDiscreteFunctionSpace<DomainSpace1Type,DomainSpace2Type,DomainSpace3Type> CombinedDiscreteFunctionSpaceType;
  typedef ISTLBlockVectorDiscreteFunction<CombinedDiscreteFunctionSpaceType> DomainFunctionType;
  typedef DomainFunctionType RangeFunctionType;
  typedef typename DomainFunctionType::RangeFieldType DomainFieldType;
  typedef typename RangeFunctionType::RangeFieldType RangeFieldType;

  typedef typename DomainFunctionType::DofStorageType domain_type;
  typedef typename RangeFunctionType::DofStorageType range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};

  typedef typename Oper11Type::RangeFunctionType RangeFunction1Type;
  typedef typename Oper21Type::RangeFunctionType RangeFunction2Type;
  typedef typename Oper31Type::RangeFunctionType RangeFunction3Type;

  // constructor
  ExtendedOperatorWrapper(const Oper11Type& oper11,const Oper12Type& oper12,const Oper21Type& oper21,const Oper13Type& oper13,
                                 const Oper31Type& oper31):
    oper11_(oper11),oper12_(oper12),oper21_(oper21),oper13_(oper13),oper31_(oper31)
  {}

  void apply(const domain_type& x,range_type& b) const
  {
    // split x into x1, x2 and x3
    DomainFunction1Type x1("x1",oper11_.domainSpace());
    DomainFunction2Type x2("x2",oper12_.domainSpace());
    DomainFunction3Type x3("x3",oper13_.domainSpace());
    std::size_t counter(0);
    for(auto& dof:dofs(x1))
    {
      dof=x[counter];
      ++counter;
    }
    for(auto& dof:dofs(x2))
    {
      dof=x[counter];
      ++counter;
    }
    for(auto& dof:dofs(x3))
    {
      dof=x[counter];
      ++counter;
    }

    // create b1, b1Temp, b2 and b3
    RangeFunction1Type b1("b1",oper11_.rangeSpace());
    RangeFunction1Type b1Temp("b1Temp",oper11_.rangeSpace());
    RangeFunction2Type b2("b2",oper21_.rangeSpace());
    RangeFunction3Type b3("b3",oper31_.rangeSpace());

    // b1 = A11*x1 + A12*x2 + A13*x3
    oper11_(x1,b1);
    oper12_(x2,b1Temp);
    b1+=b1Temp;
    oper13_(x3,b1Temp);
    b1+=b1Temp;

    // b2 = A21*x1
    oper21_(x1,b2);

    // b3 = A31*x3
    oper31_(x1,b3);

    // merge b1, b2 and b3 into b
    counter=0;
    for(auto& dof:dofs(b1))
    {
      b[counter]=dof;
      ++counter;
    }
    for(auto& dof:dofs(b2))
    {
      b[counter]=dof;
      ++counter;
    }
    for(auto& dof:dofs(b3))
    {
      b[counter]=dof;
      ++counter;
    }
  }

  void applyscaleadd(field_type alpha,const domain_type& x,range_type& b) const
  {
    range_type temp(b);
    apply(x,temp);
    const auto size(oper11_.domainSpace().size()+oper12_.domainSpace().size()+oper13_.domainSpace().size());
    for(auto i=decltype(size){0};i!=size;++i)
      b[i]+=alpha*temp[i];
  }

  private:
  const Oper11Type& oper11_;
  const Oper12Type& oper12_;
  const Oper21Type& oper21_;
  const Oper13Type& oper13_;
  const Oper31Type& oper31_;
};

}
}
#endif // DUNE_FEM_OPERATORWRAPPER_HH
