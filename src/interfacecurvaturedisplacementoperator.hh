#ifndef DUNE_FEM_INTERFACECURVATUREDISPLACEMENTOPERATOR_HH
#define DUNE_FEM_INTERFACECURVATUREDISPLACEMENTOPERATOR_HH

// deprecation warning
#warning ("WARNING : interfacecurvaturedisplacementoperator.hh is deprecated")

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include "normal.hh"

#include <vector>
#include <fstream>
#include <string>

namespace Dune
{
namespace Fem
{

template<typename LinearOperatorImp,typename BulkInterfaceGridMapperImp>
class InterfaceCurvatureDisplacementOperator:public Operator<typename LinearOperatorImp::DomainFunctionType,
                                                             typename LinearOperatorImp::RangeFunctionType>
{
  public:
  typedef LinearOperatorImp LinearOperatorType;
  typedef BulkInterfaceGridMapperImp BulkInterfaceGridMapperType;
  typedef typename LinearOperatorType::DomainFunctionType DomainFunctionType;
  typedef typename LinearOperatorType::RangeFunctionType RangeFunctionType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;
  typedef InterfaceCurvatureDisplacementOperator<LinearOperatorType,BulkInterfaceGridMapperType> ThisType;

  // constructor
  explicit InterfaceCurvatureDisplacementOperator(const DomainSpaceType& domainSpace,const RangeSpaceType& rangeSpace,
                                                  const BulkInterfaceGridMapperType& mapper):
    domainspace_(domainSpace),rangespace_(rangeSpace),mapper_(mapper),
    op_("interface curvature displacement operator",domainspace_,rangespace_)
  {}

  InterfaceCurvatureDisplacementOperator(const ThisType& )=delete;

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="interface_curvature_displacement_matrix.dat") const
  {
    std::ofstream ofs(filename);
    op_.matrix().print(ofs);
  }

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  const DomainSpaceType& domainSpace() const
  {
    return domainspace_;
  }

  const RangeSpaceType& rangeSpace() const
  {
    return rangespace_;
  }

  void assemble()
  {
    DiagonalAndNeighborStencil<DomainSpaceType,RangeSpaceType> stencil(domainspace_,rangespace_);
    op_.reserve(stencil);
    op_.clear();

    const auto domainLocalBlockSize(DomainSpaceType::localBlockSize);
    const auto rangeLocalBlockSize(RangeSpaceType::localBlockSize);
    typedef typename DomainSpaceType::RangeType DomainType;
    typedef typename RangeSpaceType::RangeType RangeType;
    std::vector<DomainType> phiDomain(domainspace_.blockMapper().maxNumDofs()*domainLocalBlockSize);
    std::vector<RangeType> phiRange(rangespace_.blockMapper().maxNumDofs()*rangeLocalBlockSize);

    // define normal functor and normal vector
    typedef typename DomainSpaceType::GridType::ctype ctype;
    constexpr auto worlddim(DomainSpaceType::GridType::dimensionworld);
    typedef Normal<ctype,worlddim> NormalType;
    NormalType normal;
    typename NormalType::NormalVectorType normalVector;

    // define selector
    constexpr auto dim(DomainType::dimension);
    typedef Selector<dim,typename DomainSpaceType::RangeFieldType> SelectorType;
    SelectorType selector;

    // perform a grid walkthrough and assemble the global matrix
    for(const auto& entity:domainspace_)
    {
      // compute normal
      const auto& faceLocalIdx(mapper_.faceLocalIdxInterface2Bulk(domainspace_.grid().leafIndexSet().index(entity)));
      normal(entity,normalVector,faceLocalIdx);

      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& domainBaseSet(localMatrix.domainBasisFunctionSet());
      const auto& rangeBaseSet(localMatrix.rangeBasisFunctionSet());

      CachingQuadrature<typename DomainSpaceType::GridPartType,0> quadrature(entity,2*domainspace_.order()+1);
      for(const auto& qp:quadrature)
      {
        domainBaseSet.evaluateAll(qp,phiDomain);
        rangeBaseSet.evaluateAll(qp,phiRange);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

        const auto rows(localMatrix.rows());
        const auto columns(localMatrix.columns());
        for(auto i=0;i!=rows;++i)
        {
          for(auto j=0;j!=columns;++j)
          {
            auto value(selector(phiRange,phiDomain,normalVector,i,j));
            value*=weight;
            localMatrix.add(i,j,value);
          }
        }
      }
    }
  }

  private:
  template<std::size_t dim,typename R>
  struct Selector;

  template<typename R>
  struct Selector<1,R>
  {
    template<typename P,typename C,typename N>
    R operator()(const P& range,const C& domain,const N& normal,std::size_t i,std::size_t j) const
    {
      return domain[j]*(range[i]*normal);
    }
  };

  template<typename R>
  struct Selector<DomainSpaceType::GridType::dimensionworld,R>
  {
    template<typename P,typename C,typename N>
    R operator()(const P& range,const C& domain,const N& normal,std::size_t i,std::size_t j) const
    {
      return range[i]*(domain[j]*normal);
    }
  };

  const DomainSpaceType& domainspace_;
  const RangeSpaceType& rangespace_;
  const BulkInterfaceGridMapperType& mapper_;
  LinearOperatorType op_;
};

}
}
#endif // DUNE_FEM_INTERFACECURVATUREDISPLACEMENTOPERATOR_HH
