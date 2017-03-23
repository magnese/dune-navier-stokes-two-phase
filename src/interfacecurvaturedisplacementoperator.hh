#ifndef DUNE_FEM_INTERFACECURVATUREDISPLACEMENTOPERATOR_HH
#define DUNE_FEM_INTERFACECURVATUREDISPLACEMENTOPERATOR_HH

// deprecation warning
#warning ("WARNING : interfacecurvaturedisplacementoperator.hh is deprecated")

#include <dune/fem/io/parameter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/linear/spoperator.hh>

#include <fstream>
#include <string>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DomainFunctionImp,typename RangeFunctionImp,typename CoupledMeshManagerImp,
         template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class InterfaceCurvatureDisplacementOperator:public Operator<DomainFunctionImp,RangeFunctionImp>
{
  public:
  typedef DomainFunctionImp DomainFunctionType;
  typedef RangeFunctionImp RangeFunctionType;
  typedef CoupledMeshManagerImp CoupledMeshManagerType;
  typedef LinearOperatorImp<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;
  typedef InterfaceCurvatureDisplacementOperator<DomainFunctionType,RangeFunctionType,CoupledMeshManagerType,LinearOperatorImp>
    ThisType;

  // constructor
  explicit InterfaceCurvatureDisplacementOperator(const DomainSpaceType& domainSpace,const RangeSpaceType& rangeSpace,
                                                  const CoupledMeshManagerType& meshManager):
    domainspace_(domainSpace),rangespace_(rangeSpace),meshmanger_(meshManager),
    op_("interface curvature displacement operator",domainspace_,rangespace_)
  {}

  InterfaceCurvatureDisplacementOperator(const ThisType& )=delete;

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="interface_curvature_displacement_matrix.dat",unsigned int offset=0) const
  {
    std::ofstream ofs(Parameter::getValue<std::string>("fem.prefix",".")+"/"+filename);
    op_.matrix().print(ofs,offset);
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

    constexpr std::size_t domainLocalBlockSize(DomainSpaceType::localBlockSize);
    constexpr std::size_t rangeLocalBlockSize(RangeSpaceType::localBlockSize);
    typedef typename DomainSpaceType::RangeType DomainType;
    typedef typename RangeSpaceType::RangeType RangeType;
    std::vector<DomainType> phiDomain(domainspace_.blockMapper().maxNumDofs()*domainLocalBlockSize);
    std::vector<RangeType> phiRange(rangespace_.blockMapper().maxNumDofs()*rangeLocalBlockSize);

    // define selector
    constexpr unsigned int dim(DomainType::dimension);
    typedef Selector<dim,typename DomainSpaceType::RangeFieldType> SelectorType;
    SelectorType selector;

    // perform a grid walkthrough and assemble the global matrix
    for(const auto& entity:domainspace_)
    {
      // compute normal
      const auto intersection(meshmanger_.correspondingInnerBulkIntersection(entity));
      const auto normalVector(intersection.centerUnitOuterNormal());

      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& domainBaseSet(localMatrix.domainBasisFunctionSet());
      const auto& rangeBaseSet(localMatrix.rangeBasisFunctionSet());

      const CachingQuadrature<typename DomainSpaceType::GridPartType,0> quadrature(entity,2*domainspace_.order()+1);
      for(const auto& qp:quadrature)
      {
        domainBaseSet.evaluateAll(qp,phiDomain);
        rangeBaseSet.evaluateAll(qp,phiRange);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());

        const auto rows(localMatrix.rows());
        const auto columns(localMatrix.columns());
        for(auto i=decltype(rows){0};i!=rows;++i)
        {
          for(auto j=decltype(columns){0};j!=columns;++j)
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
  template<unsigned int dim,typename R>
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
  const CoupledMeshManagerType& meshmanger_;
  LinearOperatorType op_;
};

}
}
#endif // DUNE_FEM_INTERFACECURVATUREDISPLACEMENTOPERATOR_HH
