#ifndef DUNE_FEM_EXTENDEDTUPLEDISCRETEFUNCTION_HH
#define DUNE_FEM_EXTENDEDTUPLEDISCRETEFUNCTION_HH

#include <algorithm>
#include <cmath>
#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/fem/function/tuplediscretefunction.hh>

namespace Dune
{
namespace Fem
{

template<typename... DiscreteFunctions>
class ExtendedTupleDiscreteFunction:public TupleDiscreteFunction<DiscreteFunctions...>
{
  public:
  typedef TupleDiscreteFunction<DiscreteFunctions...> BaseType;
  typedef ExtendedTupleDiscreteFunction<DiscreteFunctions...> ThisType;

  typedef typename BaseType::DofVectorType::FieldType FieldType;
  typedef FieldType field_type;

  using BaseType::dofVector_;
  using BaseType::scalarProductDofs;

  template<typename... Args>
  ExtendedTupleDiscreteFunction(Args&&... args):
    BaseType(std::forward<Args>(args)...)
  {}

  const ThisType& operator=(const FieldType& value)
  {
    Hybrid::forEach(typename BaseType::Sequence{},[&](auto i){std::fill(std::get<i>(*this).dbegin(),std::get<i>(*this).dend(),value);});
    return *this;
  }

  ThisType& operator=(const ThisType& other)
  {
    dofVector_=other.dofVector_;
    return *this;
  }

  FieldType dot(const ThisType& other) const
  {
    return scalarProductDofs(other);
  }

  FieldType two_norm() const
  {
    return std::sqrt(scalarProductDofs(*this));
  }
};

template<typename DiscreteFunctionType>
struct UninitializedLocalFunctionHelper
{
  static auto get(DiscreteFunctionType& df)
  {
    return std::tuple<typename DiscreteFunctionType::LocalFunctionType>(df.localFunction());
  }
};

template<typename... DiscreteFunctionsType>
struct UninitializedLocalFunctionHelper<ExtendedTupleDiscreteFunction<DiscreteFunctionsType...>>
{
  typedef ExtendedTupleDiscreteFunction<DiscreteFunctionsType...> DiscreteFunctionType;

  static auto get(DiscreteFunctionType& df)
  {
    return get(df,typename DiscreteFunctionType::Sequence{});
  }

  template<std::size_t... i>
  static auto get(DiscreteFunctionType& df,std::index_sequence<i...> )
  {
    return std::tuple<typename DiscreteFunctionsType::LocalFunctionType...>(df.template subDiscreteFunction<i>().localFunction()...);
  }
};

template<typename DiscreteFunctionType>
auto getUninitializedLocalFunctions(DiscreteFunctionType& df)
{
  return UninitializedLocalFunctionHelper<DiscreteFunctionType>::get(df);
}

}
}

#endif // DUNE_FEM_EXTENDEDTUPLEDISCRETEFUNCTION_HH
