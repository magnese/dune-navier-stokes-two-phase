#ifndef DUNE_FEM_OPERATORGLUER_HH
#define DUNE_FEM_OPERATORGLUER_HH

#include <dune/fem/io/io.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <algorithm>
#include <fstream>
#include <memory>
#include <string>
#include <tuple>

#include "extendedtuplediscretefunction.hh"

namespace Dune
{
namespace Fem
{

template<typename Op11T,typename Op12T,typename Op21T,typename Op22T>
class OperatorGluer:public Operator<ExtendedTupleDiscreteFunction<typename Op11T::DomainFunctionType,typename Op12T::DomainFunctionType>,
                                    ExtendedTupleDiscreteFunction<typename Op11T::DomainFunctionType,typename Op12T::DomainFunctionType>>
{
  private:
  typedef Op11T Op11Type;
  typedef Op12T Op12Type;
  typedef Op21T Op21Type;
  typedef Op22T Op22Type;

  typedef typename Op11Type::DomainFunctionType DomainFunction1Type;
  typedef typename Op22Type::DomainFunctionType DomainFunction2Type;

  typedef typename Op11Type::MatrixType Matrix11Type;
  typedef typename Op12Type::MatrixType Matrix12Type;
  typedef typename Op21Type::MatrixType Matrix21Type;
  typedef typename Op22Type::MatrixType Matrix22Type;

  static constexpr auto defaultCol=Matrix11Type::defaultCol;

  public:
  typedef ExtendedTupleDiscreteFunction<DomainFunction1Type,DomainFunction2Type> DiscreteFunctionType;
  typedef DiscreteFunctionType DomainFunctionType;
  typedef DomainFunctionType RangeFunctionType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  typedef DiscreteSpaceType DomainSpaceType;
  typedef DomainSpaceType RangeSpaceType;
  typedef SparseRowLinearOperator<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef OperatorGluer<Op11Type,Op12Type,Op21Type,Op22Type> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  OperatorGluer(const Op11Type& op11,const Op12Type& op12,const Op21Type& op21,const Op22Type& op22):
    mat11_(op11.systemMatrix().exportMatrix()),mat12_(op12.systemMatrix().exportMatrix()),mat21_(op21.systemMatrix().exportMatrix()),
    mat22_(op22.systemMatrix().exportMatrix()),space_(op11.domainSpace().gridPart()),op_("gluer operator",space_,space_)
  {}

  OperatorGluer(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="glued_matrix.dat",unsigned int offset=0) const
  {
    const std::string& path(Parameter::getValue<std::string>("fem.prefix","."));
    if(!directoryExists(path))
      createDirectory(path);
    std::ofstream ofs(path+"/"+filename);
    op_.exportMatrix().print(ofs,offset);
  }

  const DomainSpaceType& domainSpace() const
  {
    return space_;
  }

  const RangeSpaceType& rangeSpace() const
  {
    return space_;
  }

  void assemble()
  {
    // reserve matrix
    const auto maxNumNonZeros(std::max(mat11_.numNonZeros()+mat12_.numNonZeros(),mat21_.numNonZeros()+mat22_.numNonZeros()));
    SimpleStencil<DomainSpaceType,RangeSpaceType> stencil(maxNumNonZeros);
    op_.reserve(stencil);
    op_.clear();

    // copy matrix11 and matrix12 into operator
    const auto offset(mat11_.rows());
    auto count1(decltype(offset){0});
    auto count2(decltype(offset){0});
    auto& matrix(op_.exportMatrix());
    for(auto row=decltype(mat11_.rows()){0};row!=mat11_.rows();++row)
    {
      while(count1<(mat11_.numNonZeros()*(row+1)))
      {
        const auto entry(mat11_.realValue(count1));
        const auto col(entry.second);
        if(col!=defaultCol)
          matrix.set(row,col,entry.first);
        ++count1;
      }
      while(count2<(mat12_.numNonZeros()*(row+1)))
      {
        const auto entry(mat12_.realValue(count2));
        const auto col(entry.second);
        if(col!=defaultCol)
          matrix.set(row,col+offset,entry.first);
        ++count2;
      }
    }

    // copy matrix21 and matrix22 into operator
    count1=0;
    count2=0;
    for(auto row=decltype(mat21_.rows()){0};row!=mat21_.rows();++row)
    {
      while(count1<(mat21_.numNonZeros()*(row+1)))
      {
        const auto entry(mat21_.realValue(count1));
        const auto col(entry.second);
        if(col!=defaultCol)
          matrix.set(row+offset,col,entry.first);
        ++count1;
      }
      while(count2<(mat22_.numNonZeros()*(row+1)))
      {
        const auto entry(mat22_.realValue(count2));
        const auto col(entry.second);
        if(col!=defaultCol)
          matrix.set(row+offset,col+offset,entry.first);
        ++count2;
      }
    }
  }

  void applyDoctoring()
  {
    const auto row(mat11_.cols());
    auto& matrix(op_.exportMatrix());
    matrix.clearRow(row);
    #if USE_SYMMETRIC_DIRICHLET
    for(auto i=decltype(matrix.rows()){0};i!=matrix.rows();++i)
      matrix.set(i,row,0.0);
    #endif
    matrix.set(row,row,1.0);
  }

  private:
  const Matrix11Type& mat11_;
  const Matrix12Type& mat12_;
  const Matrix21Type& mat21_;
  const Matrix22Type& mat22_;
  const DiscreteSpaceType space_;
  LinearOperatorType op_;
};

template<typename Op11T,typename Op12T,typename Op21T,typename Op13T,typename Op31T>
class ExtendedOperatorGluer:public Operator<
  ExtendedTupleDiscreteFunction<typename Op11T::DomainFunctionType,typename Op12T::DomainFunctionType,typename Op13T::DomainFunctionType>,
  ExtendedTupleDiscreteFunction<typename Op11T::DomainFunctionType,typename Op12T::DomainFunctionType,typename Op13T::DomainFunctionType>>
{
  private:
  typedef Op11T Op11Type;
  typedef Op12T Op12Type;
  typedef Op21T Op21Type;
  typedef Op13T Op13Type;
  typedef Op31T Op31Type;

  typedef typename Op11Type::DomainFunctionType DomainFunction1Type;
  typedef typename Op12Type::DomainFunctionType DomainFunction2Type;
  typedef typename Op13Type::DomainFunctionType DomainFunction3Type;

  typedef typename Op11Type::MatrixType Matrix11Type;
  typedef typename Op12Type::MatrixType Matrix12Type;
  typedef typename Op21Type::MatrixType Matrix21Type;
  typedef typename Op13Type::MatrixType Matrix13Type;
  typedef typename Op31Type::MatrixType Matrix31Type;

  static constexpr auto defaultCol=Matrix11Type::defaultCol;

  public:
  typedef ExtendedTupleDiscreteFunction<DomainFunction1Type,DomainFunction2Type,DomainFunction3Type> DiscreteFunctionType;
  typedef DiscreteFunctionType DomainFunctionType;
  typedef DomainFunctionType RangeFunctionType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  typedef DiscreteSpaceType DomainSpaceType;
  typedef DomainSpaceType RangeSpaceType;
  typedef SparseRowLinearOperator<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef ExtendedOperatorGluer<Op11Type,Op12Type,Op21Type,Op13Type,Op31Type> ThisType;
  typedef Operator<DomainFunctionType,RangeFunctionType> BaseType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  ExtendedOperatorGluer(const Op11Type& op11,const Op12Type& op12,const Op21Type& op21,const Op13Type& op13,const Op31Type& op31):
    mat11_(op11.systemMatrix().exportMatrix()),mat12_(op12.systemMatrix().exportMatrix()),mat21_(op21.systemMatrix().exportMatrix()),
    mat13_(op13.systemMatrix().exportMatrix()),mat31_(op31.systemMatrix().exportMatrix()),
    space_(std::make_tuple(std::make_unique<typename DomainFunction1Type::DiscreteFunctionSpaceType>(op11.domainSpace().gridPart()),
                           std::make_unique<typename DomainFunction2Type::DiscreteFunctionSpaceType>(op12.domainSpace().gridPart()),
                           std::make_unique<typename DomainFunction3Type::DiscreteFunctionSpaceType>(op13.domainSpace().gridPart()))),
    op_("extended gluer operator",space_,space_)
  {}

  ExtendedOperatorGluer(const ThisType& )=delete;

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="extended_glued_matrix.dat",unsigned int offset=0) const
  {
    const std::string& path(Parameter::getValue<std::string>("fem.prefix","."));
    if(!directoryExists(path))
      createDirectory(path);
    std::ofstream ofs(path+"/"+filename);
    op_.exportMatrix().print(ofs,offset);
  }

  const DomainSpaceType& domainSpace() const
  {
    return space_;
  }

  const RangeSpaceType& rangeSpace() const
  {
    return space_;
  }

  void assemble()
  {
    // reserve matrix
    const auto maxNumNonZeros(std::max({mat11_.numNonZeros()+mat12_.numNonZeros()+mat13_.numNonZeros(),
                                        mat21_.numNonZeros(),mat31_.numNonZeros()}));
    SimpleStencil<DomainSpaceType,RangeSpaceType> stencil(maxNumNonZeros);
    op_.reserve(stencil);
    op_.clear();

    // copy matrix11, matrix12 and matrix13 into operator
    const auto offset2(mat11_.cols());
    const auto offset3(offset2+mat12_.cols());
    auto count1(decltype(offset2){0});
    auto count2(decltype(offset2){0});
    auto count3(decltype(offset2){0});
    auto& matrix(op_.exportMatrix());
    for(auto row=decltype(mat11_.rows()){0};row!=mat11_.rows();++row)
    {
      while(count1<(mat11_.numNonZeros()*(row+1)))
      {
        const auto entry(mat11_.realValue(count1));
        const auto col(entry.second);
        if(col!=defaultCol)
          matrix.set(row,col,entry.first);
        ++count1;
      }
      while(count2<(mat12_.numNonZeros()*(row+1)))
      {
        const auto entry(mat12_.realValue(count2));
        const auto col(entry.second);
        if(col!=defaultCol)
          matrix.set(row,col+offset2,entry.first);
        ++count2;
      }
      while(count3<(mat13_.numNonZeros()*(row+1)))
      {
        const auto entry(mat13_.realValue(count3));
        const auto col(entry.second);
        if(col!=defaultCol)
          matrix.set(row,col+offset3,entry.first);
        ++count3;
      }

    }

    // copy matrix21 into operator
    count1=0;
    for(auto row=decltype(mat21_.rows()){0};row!=mat21_.rows();++row)
    {
      while(count1<(mat21_.numNonZeros()*(row+1)))
      {
        const auto entry(mat21_.realValue(count1));
        const auto col(entry.second);
        if(col!=defaultCol)
          matrix.set(row+offset2,col,entry.first);
        ++count1;
      }
    }

    // copy matrix31 into operator
    count1=0;
    for(auto row=decltype(mat31_.rows()){0};row!=mat31_.rows();++row)
    {
      while(count1<(mat31_.numNonZeros()*(row+1)))
      {
        const auto entry(mat31_.realValue(count1));
        const auto col(entry.second);
        if(col!=defaultCol)
          matrix.set(row+offset3,col,entry.first);
        ++count1;
      }
    }
  }

  void applyDoctoring()
  {
    auto row(mat11_.cols());
    auto& matrix(op_.exportMatrix());
    matrix.clearRow(row);
    #if USE_SYMMETRIC_DIRICHLET
    for(auto i=decltype(matrix.rows()){0};i!=matrix.rows();++i)
      matrix.set(i,row,0.0);
    #endif
    matrix.set(row,row,1.0);
    #if PRESSURE_SPACE_TYPE != 3
    row+=mat12_.cols();
    matrix.clearRow(row);
    #if USE_SYMMETRIC_DIRICHLET
    for(auto i=decltype(matrix.rows()){0};i!=matrix.rows();++i)
      matrix.set(i,row,0.0);
    #endif
    matrix.set(row,row,1.0);
    #endif
  }

  private:
  const Matrix11Type& mat11_;
  const Matrix12Type& mat12_;
  const Matrix21Type& mat21_;
  const Matrix13Type& mat13_;
  const Matrix31Type& mat31_;
  const DiscreteSpaceType space_;
  LinearOperatorType op_;
};

}
}
#endif // DUNE_FEM_OPERATORGLUER_HH
