#ifndef DUNE_FEM_OPERATORGLUER_HH
#define DUNE_FEM_OPERATORGLUER_HH

#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <algorithm>
#include <fstream>
#include <string>

namespace Dune
{
namespace Fem
{

template<typename Op11T,typename Op12T,typename Op21T,typename Op22T>
class OperatorGluer:public Operator<
              ISTLBlockVectorDiscreteFunction<TupleDiscreteFunctionSpace<typename Op11T::DomainSpaceType,typename Op12T::DomainSpaceType>>,
              ISTLBlockVectorDiscreteFunction<TupleDiscreteFunctionSpace<typename Op11T::DomainSpaceType,typename Op12T::DomainSpaceType>>>
{
  private:
  typedef Op11T Op11Type;
  typedef Op12T Op12Type;
  typedef Op21T Op21Type;
  typedef Op22T Op22Type;

  typedef typename Op11Type::DomainFunctionType DomainFunction1Type;
  typedef typename DomainFunction1Type::DiscreteFunctionSpaceType DomainSpace1Type;
  typedef typename Op22Type::DomainFunctionType DomainFunction2Type;
  typedef typename DomainFunction2Type::DiscreteFunctionSpaceType DomainSpace2Type;

  typedef typename Op11Type::MatrixType Matrix11Type;
  typedef typename Op12Type::MatrixType Matrix12Type;
  typedef typename Op21Type::MatrixType Matrix21Type;
  typedef typename Op22Type::MatrixType Matrix22Type;

  public:
  typedef TupleDiscreteFunctionSpace<DomainSpace1Type,DomainSpace2Type> DomainSpaceType;
  typedef DomainSpaceType RangeSpaceType;
  typedef DomainSpaceType DiscreteSpaceType;
  typedef ISTLBlockVectorDiscreteFunction<DomainSpaceType> DomainFunctionType;
  typedef DomainFunctionType RangeFunctionType;
  typedef DomainFunctionType DiscreteFunctionType;
  typedef SparseRowLinearOperator<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef OperatorGluer<Op11Type,Op12Type,Op21Type,Op22Type> ThisType;
  typedef Operator<DomainSpace1Type,DomainSpace2Type> BaseType;
  typedef typename LinearOperatorType::MatrixType MatrixType;

  // constructor
  inline OperatorGluer(const Op11Type& op11,const Op12Type& op12,const Op21Type& op21,const Op22Type& op22):
    mat11_(op11.systemMatrix().matrix()),mat12_(op12.systemMatrix().matrix()),mat21_(op21.systemMatrix().matrix()),
    mat22_(op22.systemMatrix().matrix()),space_(op11.domainSpace().gridPart()),op_("gluer operator",space_,space_)
  {}

  OperatorGluer(const ThisType& )=delete;

  inline const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  // apply the operator
  virtual inline void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  // dump system matrix into file
  void print(const std::string& filename="glued_matrix.dat") const
  {
    std::ofstream ofs(filename);
    const auto rows(op_.matrix().rows());
    auto count(decltype(rows){0});
    for(auto row=decltype(rows){0};row!=rows;++row)
    {
      while(count<(op_.matrix().numNonZeros()*(row+1)))
      {
        const auto entry(op_.matrix().realValue(count));
        const auto value(entry.first);
        const auto col(entry.second);
        if((std::abs(value)>1.e-13)&&(col>-1))
          ofs<<row+1<<" "<<col+1<<" "<<value<<std::endl;
        ++count;
      }
    }
  }

  inline const DomainSpaceType& domainSpace() const
  {
    return space_;
  }

  inline const RangeSpaceType& rangeSpace() const
  {
    return space_;
  }

  void assemble() const
  {
    // reserve matrix
    const auto maxNumNonZeros(std::max(std::max(mat11_.numNonZeros(),mat12_.numNonZeros()),
                                       std::max(mat21_.numNonZeros(),mat22_.numNonZeros())));
    SimpleStencil<DomainSpaceType,RangeSpaceType> stencil(maxNumNonZeros);
    op_.reserve(stencil);
    op_.clear();

    // copy matrix11 and matrix12 into operator
    const auto offset(mat11_.rows());
    auto count1(decltype(offset){0});
    auto count2(decltype(offset){0});
    for(auto row=decltype(mat11_.rows()){0};row!=mat11_.rows();++row)
    {
      while(count1<(mat11_.numNonZeros()*(row+1)))
      {
        const auto entry(mat11_.realValue(count1));
        const auto col(entry.second);
        if(col>-1)
          op_.matrix().set(row,col,entry.first);
        ++count1;
      }
      while(count2<(mat12_.numNonZeros()*(row+1)))
      {
        const auto entry(mat12_.realValue(count2));
        const auto col(entry.second);
        if(col>-1)
          op_.matrix().set(row,col+offset,entry.first);
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
        if(col>-1)
          op_.matrix().set(row+offset,col,entry.first);
        ++count1;
      }
      while(count2<(mat22_.numNonZeros()*(row+1)))
      {
        const auto entry(mat22_.realValue(count2));
        const auto col(entry.second);
        if(col>-1)
          op_.matrix().set(row+offset,col+offset,entry.first);
        ++count2;
      }
    }
  }

  void applyDoctoring() const
  {
    const auto row(mat11_.cols());
    op_.matrix().clearRow(row);
    op_.matrix().set(row,row,1.0);
  }

  private:
  const Matrix11Type& mat11_;
  const Matrix12Type& mat12_;
  const Matrix21Type& mat21_;
  const Matrix22Type& mat22_;
  const DiscreteSpaceType space_;
  mutable LinearOperatorType op_;
};

}
}
#endif // DUNE_FEM_OPERATORGLUER_HH
