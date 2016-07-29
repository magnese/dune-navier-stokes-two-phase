#ifndef DUNE_FEM_OPERATORWRAPPER_HH
#define DUNE_FEM_OPERATORWRAPPER_HH

#include <dune/istl/operators.hh>

#include "extendedtuplediscretefunction.hh"

namespace Dune
{
namespace Fem
{

template<typename Oper11T,typename Oper12T,typename Oper21T,typename Oper22T>
class OperatorWrapper:public Dune::LinearOperator<
  ExtendedTupleDiscreteFunction<typename Oper11T::DomainFunctionType,typename Oper12T::DomainFunctionType>,
  ExtendedTupleDiscreteFunction<typename Oper11T::RangeFunctionType,typename Oper21T::RangeFunctionType>>
{
  public:
  typedef Oper11T Oper11Type;
  typedef Oper12T Oper12Type;
  typedef Oper21T Oper21Type;
  typedef Oper22T Oper22Type;

  typedef ExtendedTupleDiscreteFunction<typename Oper11T::DomainFunctionType,typename Oper12T::DomainFunctionType> domain_type;
  typedef ExtendedTupleDiscreteFunction<typename Oper11T::RangeFunctionType,typename Oper21T::RangeFunctionType> range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};

  OperatorWrapper(const Oper11Type& oper11,const Oper12Type& oper12,const Oper21Type& oper21,const Oper22Type& oper22):
    oper11_(oper11),oper12_(oper12),oper21_(oper21),oper22_(oper22)
  {}

  void apply(const domain_type& x,range_type& b) const
  {
    // create temp df
    range_type temp(b);

    // get reference sub-df's
    const auto& x1(x.template subDiscreteFunction<0>());
    const auto& x2(x.template subDiscreteFunction<1>());
    auto& b1(b.template subDiscreteFunction<0>());
    auto& b2(b.template subDiscreteFunction<1>());
    auto& temp1(temp.template subDiscreteFunction<0>());
    auto& temp2(temp.template subDiscreteFunction<1>());

    // b1 = A11*x1
    oper11_(x1,b1);

    // temp1 = A12*x2
    oper12_(x2,temp1);

    // b2 = A21*x1
    oper21_(x1,b2);

    // temp2 = A22*x2
    oper22_(x2,temp2);

    // b += temp
    b+=temp;
  }

  void applyscaleadd(field_type alpha,const domain_type& x,range_type& b) const
  {
    range_type temp(b);
    apply(x,temp);
    b.axpy(alpha,temp);
  }

  private:
  const Oper11Type& oper11_;
  const Oper12Type& oper12_;
  const Oper21Type& oper21_;
  const Oper22Type& oper22_;
};

template<typename Oper11T,typename Oper12T,typename Oper21T,typename Oper13T,typename Oper31T>
class ExtendedOperatorWrapper:public Dune::LinearOperator<
  ExtendedTupleDiscreteFunction<
    typename Oper11T::DomainFunctionType,typename Oper12T::DomainFunctionType,typename Oper13T::DomainFunctionType>,
  ExtendedTupleDiscreteFunction<
    typename Oper11T::RangeFunctionType,typename Oper21T::RangeFunctionType,typename Oper31T::RangeFunctionType>>
{
  public:
  typedef Oper11T Oper11Type;
  typedef Oper12T Oper12Type;
  typedef Oper21T Oper21Type;
  typedef Oper13T Oper13Type;
  typedef Oper31T Oper31Type;

  typedef ExtendedTupleDiscreteFunction<typename Oper11T::DomainFunctionType,typename Oper12T::DomainFunctionType,
    typename Oper13T::DomainFunctionType> domain_type;
  typedef ExtendedTupleDiscreteFunction<typename Oper11T::RangeFunctionType,typename Oper21T::RangeFunctionType,
    typename Oper31T::RangeFunctionType> range_type;
  typedef typename domain_type::field_type field_type;
  enum {category=SolverCategory::sequential};

  ExtendedOperatorWrapper(const Oper11Type& oper11,const Oper12Type& oper12,const Oper21Type& oper21,const Oper13Type& oper13,
                          const Oper31Type& oper31):
    oper11_(oper11),oper12_(oper12),oper21_(oper21),oper13_(oper13),oper31_(oper31)
  {}

  void apply(const domain_type& x,range_type& b) const
  {
    // get reference sub-df's
    const auto& x1(x.template subDiscreteFunction<0>());
    const auto& x2(x.template subDiscreteFunction<1>());
    const auto& x3(x.template subDiscreteFunction<2>());
    auto& b1(b.template subDiscreteFunction<0>());
    auto& b2(b.template subDiscreteFunction<1>());
    auto& b3(b.template subDiscreteFunction<2>());

    // b1 = A11*x1 + A12*x2 + A13*x3
    oper11_(x1,b1);
    auto temp(b1);
    oper12_(x2,temp);
    b1+=temp;
    oper13_(x3,temp);
    b1+=temp;

    // b2 = A21*x1
    oper21_(x1,b2);

    // b3 = A31*x1
    oper31_(x1,b3);
  }

  void applyscaleadd(field_type alpha,const domain_type& x,range_type& b) const
  {
    range_type temp(b);
    apply(x,temp);
    b.axpy(alpha,temp);
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
