/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef MPL_ArrayView_h
#define MPL_ArrayView_h

#include <vector>
#include "atlas/mpl/MPL.h"
#include "atlas/util/Debug.h"
#include "atlas/util/ArrayView.h"

namespace atlas {
namespace mpl {

template <typename DATA_TYPE>
class MPL_ArrayView : public ArrayView<DATA_TYPE>
{
public:
  typedef typename ArrayView<DATA_TYPE>::value_t       value_t;
  typedef typename ArrayView<DATA_TYPE>::const_value_t const_value_t;

public:
  MPL_ArrayView();

  MPL_ArrayView( const DATA_TYPE* data, const int strides[], const int shape[], const int rank,
                 const int mpl_idxpos[], const int mpl_rank );

  MPL_ArrayView( const DATA_TYPE* data, const int strides[], const int shape[], const int rank,
                 const int mpl_idxpos );

  template <int R>
    MPL_ArrayView( const ArrayView<MPL_ArrayView::value_t,R>& arrview,
                   const int mpl_idxpos[], const int mpl_rank );

  template <int R>
     MPL_ArrayView( const ArrayView<const_value_t,R>& arrview,
                    const int mpl_idxpos[], const int mpl_rank );

  int mpl_rank() const { return mpl_rank_; }

  const std::vector<int>& mpl_shape() const { return mpl_shape_; }

  const std::vector<int>& mpl_strides() const { return mpl_strides_; }

  int mpl_shape(int idx) const { return mpl_shape_[idx]; }

  int mpl_stride(int idx) const { return mpl_strides_[idx]; }

  int var_rank() const { return var_rank_; }

  const std::vector<int>& var_shape() const { return var_shape_; }

  const std::vector<int>& var_strides() const { return var_strides_; }

  int var_shape(int idx) const { return var_shape_[idx]; }

  int var_stride(int idx) const { return var_strides_[idx]; }

private:
  void constructor(const int mpl_idxpos[], const int mpl_rank);

private:
  std::vector<int> mpl_strides_;
  std::vector<int> mpl_shape_;
  int mpl_rank_;
  int mpl_size_;
  std::vector<int> var_strides_;
  std::vector<int> var_shape_;
  int var_rank_;
  int var_size_;
};

///////////////////////////////////////////////////////////////////////////

template <typename DATA_TYPE>
MPL_ArrayView<DATA_TYPE>::MPL_ArrayView() : ArrayView<DATA_TYPE>() {}

template <typename DATA_TYPE>
MPL_ArrayView<DATA_TYPE>::MPL_ArrayView( const DATA_TYPE* data, const int strides[], const int shape[], const int rank,
               const int mpl_idxpos[], const int mpl_rank ) :
  ArrayView<DATA_TYPE>(data,strides,shape,rank)
{
  constructor(mpl_idxpos,mpl_rank);
}

template <typename DATA_TYPE>
MPL_ArrayView<DATA_TYPE>::MPL_ArrayView( const DATA_TYPE* data, const int strides[], const int shape[], const int rank,
               const int mpl_idxpos ) :
  ArrayView<DATA_TYPE>(data,strides,shape,rank)
{
  constructor(&mpl_idxpos,1);
}

template <typename DATA_TYPE>
template <int R>
MPL_ArrayView<DATA_TYPE>::MPL_ArrayView( const ArrayView<value_t,R>& arrview,
               const int mpl_idxpos[], const int mpl_rank ) :
  ArrayView<DATA_TYPE>(arrview.data(),arrview.strides(),arrview.shape(),arrview.rank())
{
  constructor(mpl_idxpos,mpl_rank);
}

template <typename DATA_TYPE>
template <int R>
MPL_ArrayView<DATA_TYPE>::MPL_ArrayView( const ArrayView<const_value_t,R>& arrview,
               const int mpl_idxpos[], const int mpl_rank ) :
  ArrayView<DATA_TYPE>(arrview.data(),arrview.strides(),arrview.shape(),arrview.rank())
{
  constructor(mpl_idxpos,mpl_rank);
}

template <typename DATA_TYPE>
void MPL_ArrayView<DATA_TYPE>::constructor(const int mpl_idxpos[], const int mpl_rank)
{
  if( ArrayView<DATA_TYPE>::rank() == 1 )
  {
    mpl_rank_ = 1;
    mpl_strides_.assign(1,ArrayView<DATA_TYPE>::stride(0));
    mpl_shape_.assign(1,ArrayView<DATA_TYPE>::shape(0));
    var_rank_ = 1;
    var_strides_.assign(1,1);
    var_shape_.assign(1,1);
  }
  else
  {
    mpl_rank_ = mpl_rank;
    var_rank_ = ArrayView<DATA_TYPE>::rank()-mpl_rank_;
    mpl_strides_.reserve(mpl_rank_);
    mpl_shape_.reserve(mpl_rank_);
    var_strides_.reserve(var_rank_);
    var_shape_.reserve(var_rank_);
    mpl_size_ = 0;
    var_size_ = 0;
    int nvar(0);
    int nmpl(0);
    for( int jpos=0; jpos<ArrayView<DATA_TYPE>::rank(); ++jpos )
    {
      bool isvar(true);
      for( int jmpl=0; jmpl<mpl_rank_; ++jmpl )
      {
        if( jpos == mpl_idxpos[jmpl] )
        {
          isvar=false;
          break;
        }
      }
      if( isvar )
      {
        var_strides_.push_back(ArrayView<DATA_TYPE>::stride(jpos));
        var_shape_.push_back(ArrayView<DATA_TYPE>::shape(jpos));
        var_size_ += var_shape_[nvar];
        ++nvar;
      }
      else
      {
        mpl_strides_.push_back(ArrayView<DATA_TYPE>::stride(jpos));
        mpl_shape_.push_back(ArrayView<DATA_TYPE>::shape(jpos));
        mpl_size_ += mpl_shape_[nmpl];
        ++nmpl;
      }
    }
  }
}

} // namespace mpl
} // namespace atlas

#endif // MPL_ArrayView_h
