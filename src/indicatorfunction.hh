#ifndef DUNE_FEM_INDICATORFUNCTION_HH
#define DUNE_FEM_INDICATORFUNCTION_HH

#include <vector>

#include <dune/fem/gridpart/filter/threadfilter.hh>

namespace Dune
{
namespace Fem
{

template<typename GridPartImp>
class InnerBulkGridFilter:public ThreadFilter<GridPartImp,std::vector<int>>
{
  public:
  typedef GridPartImp GridPartType;
  typedef ThreadFilter<GridPartType,std::vector<int>> BaseType;
  typedef InnerBulkGridFilter<GridPartType> ThisType;

  InnerBulkGridFilter(const GridPartType& gridPart,const std::vector<int>& elementIDs):
    BaseType(gridPart,elementIDs,1)
  {}
};

template<typename GridPartImp>
class OuterBulkGridFilter:public ThreadFilter<GridPartImp,std::vector<int>>
{
  public:
  typedef GridPartImp GridPartType;
  typedef ThreadFilter<GridPartType,std::vector<int>> BaseType;
  typedef InnerBulkGridFilter<GridPartType> ThisType;

  OuterBulkGridFilter(const GridPartType& gridPart,const std::vector<int>& elementIDs):
    BaseType(gridPart,elementIDs,2)
  {}
};

}
}

#endif // DUNE_FEM_INDICATORFUNCTION_HH
