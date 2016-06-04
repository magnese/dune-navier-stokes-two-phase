#ifndef DUNE_FEM_EXTENDEDTUPLEDISCRETEFUNCTION_HH
#define DUNE_FEM_EXTENDEDTUPLEDISCRETEFUNCTION_HH

#include <algorithm>
#include <cmath>
#include <utility>

#include <dune/fem/common/tupleforeach.hh>
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
    for_each(*this,[&value](auto& df,auto ){std::fill(df.dbegin(),df.dend(),value);});
    return *this;
  }

  ThisType& operator=(const ThisType& other)
  {
    dofVector_ = other.dofVector_;
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

}
}

#endif // DUNE_FEM_EXTENDEDTUPLEDISCRETEFUNCTION_HH
