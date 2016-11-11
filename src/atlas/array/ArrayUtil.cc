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

ArraySpec::ArraySpec( const ArrayShape& shape )
{
  rank_ = shape_.size();
  size_ = 1;
  shape_.resize(rank_);
  strides_.resize(rank_);
  for( int j=rank_-1; j>=0; --j ) {
    shape_[j] = shape[j];
    strides_[j] = size_;
    size_ *= shape_[j];
  }
  contiguous_ = true;
};

ArraySpec::ArraySpec( const ArrayShape& _shape, const ArrayStrides& _strides )
{
  if( _shape.size() != _strides.size() )
    throw eckit::BadParameter("dimensions of shape and stride don't match", Here());

  shape_=_shape;
  strides_=_strides;
  rank_ = shape_.size();
  size_ = 1;
  for( size_t n=0; n<rank_; ++n )
  {
    size_ *= shape_[n];
  }
  contiguous_ = (size_ == shape_[0]*strides_[0] ? true : false);
}

const std::vector<int>& ArraySpec::shapef() const
{
  if( shapef_.empty() )
  {
    shapef_.resize(shape().size());
    std::reverse_copy( shape().begin(), shape().end(), shapef_.begin() );
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
