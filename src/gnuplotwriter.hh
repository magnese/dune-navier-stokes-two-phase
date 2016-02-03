#ifndef DUNE_FEM_GNUPLOTWRITER_HH
#define DUNE_FEM_GNUPLOTWRITER_HH

#include <fstream>
#include <list>
#include <string>
#include <utility>

#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

// generic gnuplot writer
struct GnuplotWriter
{
  GnuplotWriter(std::string&& fileName):
    filename_("/"+fileName+".dat")
  {}

  std::list<std::pair<double,double>>& get()
  {
    return values_;
  }

  const std::list<std::pair<double,double>>& get() const
  {
    return values_;
  }

  void add(double&& first,double&& second)
  {
    if(values_.size()==0)
      filename_=Parameter::getValue<std::string>("fem.prefix",".")+filename_;
    values_.emplace_back(first,second);
  }

  ~GnuplotWriter()
  {
    if(values_.size()!=0)
    {
      std::ofstream file(filename_);
      for(const auto& value:values_)
        file<<value.first<<" "<<value.second<<std::endl;
    }
  }

  std::string filename_;
  std::list<std::pair<double,double>> values_;
};

}
}

#endif // DUNE_FEM_GNUPLOTWRITER_HH
