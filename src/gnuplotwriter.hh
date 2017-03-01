#ifndef DUNE_FEM_GNUPLOTWRITER_HH
#define DUNE_FEM_GNUPLOTWRITER_HH

#include <fstream>
#include <iomanip>
#include <list>
#include <string>
#include <tuple>
#include <utility>

#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

// generic gnuplot writer
struct GnuplotWriter
{
  typedef std::list<std::tuple<double,double,bool>> ListType;

  GnuplotWriter(const std::string& fileName,unsigned int precision=6):
    filename_(Parameter::getValue<std::string>("fem.prefix",".")+"/"+fileName+".dat"),precision_(precision)
  {}

  void add(double first,double second,bool leaveEmptyRow=false)
  {
    values_.emplace_back(first,second,leaveEmptyRow);
  }

  bool isEmpty() const
  {
    return values_.size()==0;
  }

  void normalize(double coeff)
  {
    for(auto& value:values_)
      std::get<1>(value)/=coeff;
  }

  std::pair<double,double> firstValue() const
  {
    return {std::get<0>(values_.front()),std::get<1>(values_.front())};
  }

  std::pair<double,double> lastValue() const
  {
    return {std::get<0>(values_.back()),std::get<1>(values_.back())};
  }

  std::pair<double,double> minValue() const
  {
    auto minVal(firstValue());
    for(const auto& value:values_)
      if(std::get<1>(value)<minVal.second)
        minVal={std::get<0>(value),std::get<1>(value)};
    return minVal;
  }

  std::pair<double,double> maxValue() const
  {
    auto maxVal(firstValue());
    for(const auto& value:values_)
      if(std::get<1>(value)>maxVal.second)
        maxVal={std::get<0>(value),std::get<1>(value)};
    return maxVal;
  }

  ~GnuplotWriter()
  {
    finalize();
  }

  void finalize() const
  {
    if(!isEmpty())
    {
      std::ofstream file(filename_);
      file<<std::setprecision(precision_);
      for(const auto& value:values_)
      {
        file<<std::get<0>(value)<<" "<<std::get<1>(value)<<"\n";
        if(std::get<2>(value))
          file<<"\n";
      }
    }
  }

  std::string filename_;
  unsigned int precision_;
  ListType values_;
};

}
}

#endif // DUNE_FEM_GNUPLOTWRITER_HH
