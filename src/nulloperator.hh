#ifndef DUNE_FEM_NULLOPERATOR_HH
#define DUNE_FEM_NULLOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/linear/spoperator.hh>

#include <string>
#include <type_traits>

namespace Dune
{
namespace Fem
{

template<typename DomainFunctionImp,typename RangeFunctionImp,
         template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class NullOperator:public Operator<DomainFunctionImp,RangeFunctionImp>
{
  public:
  typedef DomainFunctionImp DomainFunctionType;
  typedef RangeFunctionImp RangeFunctionType;
  typedef LinearOperatorImp<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef NullOperator<DomainFunctionType,RangeFunctionType,LinearOperatorImp> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  explicit NullOperator(const DomainSpaceType& domainSpace,const RangeSpaceType& rangeSpace):
    domainspace_(domainSpace),rangespace_(rangeSpace),op_("null operator",domainspace_,rangespace_)
  {}

  NullOperator(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  virtual void operator()(const DomainFunctionType& ,RangeFunctionType& w) const
  {
    w.clear();
  }

  template<typename... Args>
  void print(Args&&... ) const
  {}

  const DomainSpaceType& domainSpace() const
  {
    return domainspace_;
  }

  const RangeSpaceType& rangeSpace() const
  {
    return rangespace_;
  }

  void assemble()
  {}

  private:
  const DomainSpaceType& domainspace_;
  const RangeSpaceType& rangespace_;
  LinearOperatorType op_;
};


template<typename >
struct isNullOperator:std::false_type
{};

template<typename... Args>
struct isNullOperator<NullOperator<Args...>>:std::true_type
{};

}
}

#endif // DUNE_FEM_NULLOPERATOR_HH
