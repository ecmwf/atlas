/*
 * (C) Copyright 1996-2015 ECMWF.
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
#include "atlas/mpi/mpi.h"
#include "atlas/util/Debug.h"
#include "atlas/util/ArrayView.h"

namespace atlas {
namespace mpl {

template <typename DATA_TYPE>
class MPL_ArrayView : public ArrayView<DATA_TYPE>
{
public:
  typedef typename ArrayView<DATA_TYPE>::value_type       value_type;

public:
  MPL_ArrayView();

  MPL_ArrayView( const DATA_TYPE* data, const size_t strides[], const size_t shape[], const size_t rank,
                 const size_t mpl_idxpos[], const size_t mpl_rank );

  MPL_ArrayView( const DATA_TYPE* data, const size_t strides[], const size_t shape[], const size_t rank,
                 const size_t mpl_idxpos );

  template <int R>
    MPL_ArrayView( const ArrayView<MPL_ArrayView::value_type,R>& arrview,
                   const size_t mpl_idxpos[], const size_t mpl_rank );

  template <int R>
     MPL_ArrayView( const ArrayView<const value_type,R>& arrview,
                    const size_t mpl_idxpos[], const size_t mpl_rank );

  int mpl_rank() const { return mpl_rank_; }

  const std::vector<size_t>& mpl_shape() const { return mpl_shape_; }

  const std::vector<size_t>& mpl_strides() const { return mpl_strides_; }

  size_t mpl_shape(size_t idx) const { return mpl_shape_[idx]; }

  size_t mpl_stride(size_t idx) const { return mpl_strides_[idx]; }

  size_t var_rank() const { return var_rank_; }

  const std::vector<size_t>& var_shape() const { return var_shape_; }

  const std::vector<size_t>& var_strides() const { return var_strides_; }

  size_t var_shape(size_t idx) const { return var_shape_[idx]; }

  size_t var_stride(size_t idx) const { return var_strides_[idx]; }

private:
  void constructor(const size_t mpl_idxpos[], const size_t mpl_rank);

private:
  std::vector<size_t> mpl_strides_;
  std::vector<size_t> mpl_shape_;
  size_t mpl_rank_;
  size_t mpl_size_;
  std::vector<size_t> var_strides_;
  std::vector<size_t> var_shape_;
  size_t var_rank_;
  size_t var_size_;
};

///////////////////////////////////////////////////////////////////////////

template <typename DATA_TYPE>
MPL_ArrayView<DATA_TYPE>::MPL_ArrayView() : ArrayView<DATA_TYPE>() {}

template <typename DATA_TYPE>
MPL_ArrayView<DATA_TYPE>::MPL_ArrayView( const DATA_TYPE* data, const size_t strides[], const size_t shape[], const size_t rank,
               const size_t mpl_idxpos[], const size_t mpl_rank ) :
  ArrayView<DATA_TYPE>(data,rank,shape,strides)
{
  constructor(mpl_idxpos,mpl_rank);
}

template <typename DATA_TYPE>
MPL_ArrayView<DATA_TYPE>::MPL_ArrayView( const DATA_TYPE* data, const size_t strides[], const size_t shape[], const size_t rank,
               const size_t mpl_idxpos ) :
  ArrayView<DATA_TYPE>(data,rank,shape,strides)
{
  constructor(&mpl_idxpos,1);
}

template <typename DATA_TYPE>
template <int R>
MPL_ArrayView<DATA_TYPE>::MPL_ArrayView( const ArrayView<value_type,R>& arrview,
               const size_t mpl_idxpos[], const size_t mpl_rank ) :
  ArrayView<DATA_TYPE>(arrview.data(),arrview.rank(),arrview.shape(),arrview.strides())
{
  constructor(mpl_idxpos,mpl_rank);
}

template <typename DATA_TYPE>
template <int R>
MPL_ArrayView<DATA_TYPE>::MPL_ArrayView( const ArrayView<const value_type,R>& arrview,
               const size_t mpl_idxpos[], const size_t mpl_rank ) :
  ArrayView<DATA_TYPE>(arrview.data(),arrview.rank(),arrview.shape(),arrview.strides())
{
  constructor(mpl_idxpos,mpl_rank);
}

template <typename DATA_TYPE>
void MPL_ArrayView<DATA_TYPE>::constructor(const size_t mpl_idxpos[], const size_t mpl_rank)
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
    size_t nvar(0);
    size_t nmpl(0);
    for( size_t jpos=0; jpos<ArrayView<DATA_TYPE>::rank(); ++jpos )
    {
      bool isvar(true);
      for( size_t jmpl=0; jmpl<mpl_rank_; ++jmpl )
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
