#ifndef DUNE_FEM_GNUPLOTWRITER_HH
#define DUNE_FEM_GNUPLOTWRITER_HH

#include <fstream>
#include <iomanip>
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
  GnuplotWriter(const std::string& fileName,unsigned int precision=6):
    filename_("/"+fileName+".dat"),filenameset_(false),precision_(precision)
  {}

  void setPrecision(unsigned int precision)
  {
    precision_=precision;
  }

  void setFileName(const std::string& fileName)
  {
    filename_=Parameter::getValue<std::string>("fem.prefix",".")+"/"+fileName+".dat";
    filenameset_=true;
  }

  std::list<std::pair<double,double>>& get()
  {
    return values_;
  }

  const std::list<std::pair<double,double>>& get() const
  {
    return values_;
  }

  void add(double first,double second,bool leaveEmptyRow=false)
  {
    if(!filenameset_)
    {
      filename_=Parameter::getValue<std::string>("fem.prefix",".")+filename_;
      filenameset_=true;
    }
    values_.emplace_back(first,second);
    returns_.emplace_back(leaveEmptyRow);
  }

  ~GnuplotWriter()
  {
    finalize();
  }

  void finalize()
  {
    if(values_.size()!=0)
    {
      std::ofstream file(filename_);
      file<<std::setprecision(precision_);
      auto returnsIt(returns_.begin());
      for(const auto& value:values_)
      {
        file<<value.first<<" "<<value.second<<std::endl;
        if(*returnsIt)
          file<<std::endl;
        ++returnsIt;
      }
    }
  }

  std::string filename_;
  bool filenameset_;
  unsigned int precision_;
  std::list<std::pair<double,double>> values_;
  std::list<bool> returns_;
};

}
}

#endif // DUNE_FEM_GNUPLOTWRITER_HH
