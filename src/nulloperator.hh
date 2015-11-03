#ifndef DUNE_FEM_NULLOPERATOR_HH
#define DUNE_FEM_NULLOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/common/unused.hh>

#include <string>
#include <fstream>
#include <type_traits>

namespace Dune
{
namespace Fem
{

template<class LinearOperatorImp,typename... Args>
class NullOperator:public Operator<typename LinearOperatorImp::DomainFunctionType,typename LinearOperatorImp::RangeFunctionType>
{
  public:
  typedef LinearOperatorImp LinearOperatorType;
  typedef typename LinearOperatorType::DomainFunctionType DomainFunctionType;
  typedef typename LinearOperatorType::RangeFunctionType RangeFunctionType;
  typedef NullOperator<LinearOperatorType,Args...> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  template<typename... Argss>
  explicit NullOperator(const DomainSpaceType& domainSpace,const RangeSpaceType& rangeSpace,const Argss&... ):
    domainspace_(domainSpace),rangespace_(rangeSpace),op_("null operator",domainspace_,rangespace_)
  {}

  NullOperator(const ThisType& )=delete;

  inline LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  // apply the operator
  virtual inline void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    DUNE_UNUSED_PARAMETER(u);
    w.clear();
  }

  // dump system matrix into file
  inline void print(const std::string& filename="null_matrix.dat") const
  {
    std::ofstream ofs(filename);
    op_.matrix().print(ofs,1);
    ofs.close();
  }

  inline const DomainSpaceType& domainSpace() const
  {
    return domainspace_;
  }

  inline const RangeSpaceType& rangeSpace() const
  {
    return rangespace_;
  }

  template<typename... Argss>
  void assemble(const Argss&... ) const
  {
    DiagonalAndNeighborStencil<DomainSpaceType,RangeSpaceType> stencil(domainspace_,rangespace_);
    op_.reserve(stencil);
    op_.clear();
  }

  private:
  const DomainSpaceType& domainspace_;
  const RangeSpaceType& rangespace_;
  mutable LinearOperatorType op_;
};


template<typename>
struct isNullOperator:std::false_type
{};

template<typename... Args>
struct isNullOperator<NullOperator<Args...>>:std::true_type
{};

}
}

#endif // DUNE_FEM_NULLOPERATOR_HH
