#ifndef DUNE_FEM_BUBBLESTATISTICS_HH
#define DUNE_FEM_BUBBLESTATISTICS_HH

#include <cmath>
#include <iostream>
#include <string>

#include "gnuplotwriter.hh"

namespace Dune
{
namespace Fem
{

// compute various bubble statistics
class BubbleStatistics
{
  public:
  BubbleStatistics(unsigned int precision=6):
    circularitywriter_("circularity",precision),barycenterwriter_("barycenter",precision),velocitywriter_("velocity",precision)
  {}

  ~BubbleStatistics()
  {
    printInfo();
  }

  template<typename FluidStateType,typename TimeProviderType>
  void add(const FluidStateType& fluidState,const TimeProviderType& timeProvider)
  {
    // compute circularity
    const auto bulkInnerVolume(fluidState.meshManager().bulkInnerVolume());
    const auto interfaceLength(fluidState.meshManager().interfaceLength());
    constexpr auto worlddim(FluidStateType::BulkGridType::dimensionworld);
    circularitywriter_.add(timeProvider.time(),circularity<worlddim>(bulkInnerVolume,interfaceLength));
    // compute height barycenter
    #if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
    barycenterwriter_.add(timeProvider.time(),verticalComponentInnerIntegrationHostEntity(
                                                fluidState.bulkGrid().coordFunction().discreteFunction(),bulkInnerVolume,
                                                fluidState.bulkInnerGridPart()));
    #endif
    // compute average rising velocity
    velocitywriter_.add(timeProvider.time(),verticalComponentInnerIntegration(fluidState.velocity(),bulkInnerVolume,
                                                                              fluidState.bulkInnerGridPart()));
  }

  void printInfo() const
  {
    if(!circularitywriter_.isEmpty())
    {
      const auto minCircularity(circularitywriter_.minValue());
      std::cout<<"Minimum circularity = "<<minCircularity.second<<" (time = "<<minCircularity.first<<" s)."<<std::endl;
      const auto maxRisingVelocity(velocitywriter_.maxValue());
      std::cout<<"Maximum average rising velocity  = "<<maxRisingVelocity.second<<" (time = "<<maxRisingVelocity.first<<" s)."<<std::endl;
      #if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
      const auto finalBarycenter(barycenterwriter_.lastValue());
      std::cout<<"Final height barycenter  = "<<finalBarycenter.second<<" (time = "<<finalBarycenter.first<<" s)."<<std::endl;
      #endif
    }
  }

  private:
  GnuplotWriter circularitywriter_;
  GnuplotWriter barycenterwriter_;
  GnuplotWriter velocitywriter_;

  template<int worlddim>
  double circularity(double bulkInnerVolume,double interfaceLength) const
  {
    // compute (perimeter_area-equivalent_sphere)/perimeter_bubble
    constexpr double d(static_cast<double>(worlddim));
    const auto equivalentCircleLength(std::pow(M_PI,1.0/d)*std::pow(2.0*d*bulkInnerVolume,(d-1.0)/d));
    return equivalentCircleLength/interfaceLength;
  }

  template<typename DiscreteFunctionType,typename BulkInnerGridPartType>
  double verticalComponentInnerIntegration(const DiscreteFunctionType& f,double bulkInnerVolume,
                                           const BulkInnerGridPartType& bulkInnerGridPart) const
  {
    // compute \int_{\Omega_-}(\vec f * \vec e_d) / \int_{\Omega_-}( 1 )
    typedef typename DiscreteFunctionType::RangeType RangeType;
    RangeType e_d(0.0);
    e_d[DiscreteFunctionType::GridType::dimensionworld-1]=1.0;
    double integral(0.0);
    for(const auto& entity:elements(bulkInnerGridPart))
    {
      const auto fLocal(f.localFunction(entity));
      CachingQuadrature<typename DiscreteFunctionType::GridPartType,0> quadrature(entity,2*f.space().order()+1);
      for(const auto& qp:quadrature)
      {
        RangeType fValue;
        fLocal.evaluate(qp.position(),fValue);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        integral+=((fValue*e_d)*weight);
      }
    }
    return integral/bulkInnerVolume;
  }

  #if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  template<typename DiscreteFunctionType,typename BulkInnerGridPartType>
  double verticalComponentInnerIntegrationHostEntity(const DiscreteFunctionType& f,double bulkInnerVolume,
                                                     const BulkInnerGridPartType& bulkInnerGridPart) const
  {
    // compute \int_{\Omega_-}(\vec f * \vec e_d) / \int_{\Omega_-}( 1 )
    typedef typename DiscreteFunctionType::RangeType RangeType;
    RangeType e_d(0.0);
    e_d[DiscreteFunctionType::GridType::dimensionworld-1]=1.0;
    double integral(0.0);
    for(const auto& entity:elements(bulkInnerGridPart))
    {
      const auto& hostEntity(entity.impl().hostEntity());
      const auto fLocal(f.localFunction(hostEntity));
      CachingQuadrature<typename DiscreteFunctionType::GridPartType,0> quadrature(hostEntity,2*f.space().order()+1);
      for(const auto& qp:quadrature)
      {
        RangeType fValue;
        fLocal.evaluate(qp.position(),fValue);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        integral+=((fValue*e_d)*weight);
      }
    }
    return integral/bulkInnerVolume;
  }
  #endif
};

}
}

#endif // DUNE_FEM_BUBBLE_STATISTICS_HH
