/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include <vector>
#include "eckit/exception/Exceptions.h"
#include "atlas/array/ArrayUtil.h"

namespace atlas {
namespace array {

ArraySpec::ArraySpec(): 
    size_(0),
    rank_(0),
    contiguous_(true),
    default_layout_(true)
{
}

ArraySpec::ArraySpec( const ArrayShape& shape )
{
  rank_ = shape.size();
  size_ = 1;
  shape_.resize(rank_);
  strides_.resize(rank_);
  layout_.resize(rank_);
  for( int j=rank_-1; j>=0; --j ) {
    shape_[j]   = shape[j];
    strides_[j] = size_;
    layout_[j]  = j;
    size_ *= shape_[j];
  }
  contiguous_ = true;
  default_layout_ = true;
};

ArraySpec::ArraySpec( const ArrayShape& shape, const ArrayStrides& strides )
{
  if( shape.size() != strides.size() )
    throw eckit::BadParameter("dimensions of shape and stride don't match", Here());

  rank_ = shape.size();
  size_ = 1;
  shape_.resize(rank_);
  strides_.resize(rank_);
  layout_.resize(rank_);
  for( int j=rank_-1; j>=0; --j ) {
    shape_[j]   = shape[j];
    strides_[j] = strides[j];
    layout_[j]  = j;
    size_ *= shape_[j];
  }
  contiguous_ = (size_ == shape_[0]*strides_[0] ? true : false);
  default_layout_ = true;
}

ArraySpec::ArraySpec( const ArrayShape& shape, const ArrayStrides& strides, const ArrayLayout& layout )
{
  if( shape.size() != strides.size() )
    throw eckit::BadParameter("dimensions of shape and stride don't match", Here());

  rank_ = shape.size();
  size_ = 1;
  shape_.resize(rank_);
  strides_.resize(rank_);
  layout_.resize(rank_);
  default_layout_ = true;
  for( int j=rank_-1; j>=0; --j ) {
    shape_[j]   = shape[j];
    strides_[j] = strides[j];
    layout_[j]  = layout[j];
    size_ *= shape_[j];
    if( layout_[j] != j ) default_layout_ = false;
  }
  contiguous_ = (size_ == shape_[0]*strides_[0] ? true : false);
}

const std::vector<int>& ArraySpec::shapef() const
{
  if( shapef_.empty() ) {
    shapef_.resize(rank_);
    for( size_t j=0; j<rank_; ++j ) {
      shapef_[j] = shape_[rank_-1- layout_[j] ];
    }
  }
  return shapef_;
}

const std::vector<int>& ArraySpec::stridesf() const
{
  if( stridesf_.empty() )
  {
    stridesf_.resize(strides().size());
    std::reverse_copy( strides().begin(), strides().end(), stridesf_.begin() );
  }
  return stridesf_;
}


} // namespace array
} // namespace atlas
