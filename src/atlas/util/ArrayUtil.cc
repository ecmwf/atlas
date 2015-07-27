/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include <vector>
#include "atlas/util/ArrayUtil.h"

namespace atlas {

ArraySpec::ArraySpec( const ArrayShape& _shape )
{
  shape_=_shape;
  strides_.resize(shape_.size());
  strides_[shape_.size()-1] = 1;
  for( long n=shape_.size()-2; n>=0; --n )
  {
    strides_[n] = strides_[n+1]*shape_[n+1];
  }
  rank_ = shape_.size();
  size_ = shape_[0]*strides_[0];
};

const std::vector<int>& ArraySpec::shapef() const
{
  if( shapef_.empty() )
  {
    shapef_.resize(shape().size());
    std::reverse_copy( shape().begin(), shape().end(), shapef_.begin() );
  }
  return shapef_;
}

} // namespace atlas
