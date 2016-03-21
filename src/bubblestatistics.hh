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
static class BubbleStatistics
{
  public:
  BubbleStatistics():
    circularitywriter_("circularity"),barycenterwriter_("barycenter"),velocitywriter_("velocity")
  {}

  ~BubbleStatistics()
  {
    if(circularitywriter_.get().size()!=0)
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
    barycenterwriter_.add(timeProvider.time(),verticalComponentInnerIntegration(
                                                fluidState.meshManager().bulkGrid().coordFunction().discreteFunction(),bulkInnerVolume,
                                                fluidState.meshManager().bulkIndicatorFunction()));
    // compute average rising velocity
    velocitywriter_.add(timeProvider.time(),verticalComponentInnerIntegration(fluidState.velocity(),bulkInnerVolume,
                                                                              fluidState.meshManager().bulkIndicatorFunction()));
  }

  void printInfo() const
  {
    // print min circularity
    const auto& circularities(circularitywriter_.get());
    auto minCircularity(circularities.front());
    for(const auto& entry:circularities)
      if(entry.second<minCircularity.second)
        minCircularity=entry;
    std::cout<<"Minimum circularity = "<<minCircularity.second<<" (time = "<<minCircularity.first<<" s)."<<std::endl;
    // print max average rising velocity
    const auto& velocities(velocitywriter_.get());
    auto maxRisingVelocity(velocities.front());
    for(const auto& entry:velocities)
      if(entry.second>maxRisingVelocity.second)
        maxRisingVelocity=entry;
    std::cout<<"Maximum average rising velocity  = "<<maxRisingVelocity.second<<" (time = "<<maxRisingVelocity.first<<" s)."<<std::endl;
    // print final height barycenter
    const auto& finalBarycenter(barycenterwriter_.get().back());
    std::cout<<"Final height barycenter  = "<<finalBarycenter.second<<" (time = "<<finalBarycenter.first<<" s)."<<std::endl;
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

  template<typename DiscreteFunctionType,typename BulkIndicatorFunctionType>
  double verticalComponentInnerIntegration(const DiscreteFunctionType& f,double bulkInnerVolume,
                                           const BulkIndicatorFunctionType& indicator) const
  {
    // compute \int_{\Omega_-}(\vec f * \vec e_d) / \int_{\Omega_-}( 1 )
    const auto& space(f.space());
    typedef typename DiscreteFunctionType::RangeType RangeType;
    RangeType e_d(0.0);
    e_d[DiscreteFunctionType::GridType::dimensionworld-1]=1.0;
    double integral(0.0);
    for(const auto& entity:space)
      if(indicator.isInner(entity))
      {
        const auto fLocal(f.localFunction(entity));
        CachingQuadrature<typename DiscreteFunctionType::GridPartType,0> quadrature(entity,2*space.order()+1);
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

} bubbleStatistics;

}
}

#endif // DUNE_FEM_BUBBLE_STATISTICS_HH
