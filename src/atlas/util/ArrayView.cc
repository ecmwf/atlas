/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <stdexcept>

#include "atlas/Array.h"
#include "atlas/util/ArrayView.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

#define TEMPLATE_SPECIALIZATION( DATA_TYPE ) \
template<>\
ArrayView <DATA_TYPE, 0 >::ArrayView( const Array& array ) : data_( const_cast<DATA_TYPE*>(array.data<DATA_TYPE>()) ) \
{ \
  rank_ = array.strides().size(); \
  strides_ = array.strides(); \
  shape_ = array.shape(); \
  size_ = std::accumulate(shape_.data(),shape_.data()+rank_,1,std::multiplies<int>()); \
} \
template<>\
ArrayView <DATA_TYPE, 1 >::ArrayView( const Array& array ) : data_( const_cast<DATA_TYPE*>(array.data<DATA_TYPE>()) ) \
{ \
  ASSERT( array.shape().size() >= 1 ); \
  strides_[0]=array.stride(0);  shape_[0]=array.shape(0); \
  ASSERT(size() == array.size()); \
} \
template<> \
ArrayView<DATA_TYPE,2>::ArrayView( const Array& array ) : data_( const_cast<DATA_TYPE*>(array.data<DATA_TYPE>()) ) \
{ \
  ASSERT( array.rank() > 0 ); \
  if( array.rank() == 1  ) { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = 1;               strides_[1] = 1; \
  } else if ( array.rank() == 2 ) { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = array.shape(1);  strides_[1] = array.stride(1); \
  } else { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = array.stride(0); strides_[1] = 1; \
  } \
  size_ = shape_[0]*shape_[1]; \
  ASSERT(size() == array.size()); \
} \
template<> \
ArrayView<DATA_TYPE,3>::ArrayView( const Array& array ) : data_( const_cast<DATA_TYPE*>(array.data<DATA_TYPE>()) ) \
{ \
  ASSERT( array.rank() > 0 ); \
  if( array.rank() == 1  ) { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = 1;               strides_[1] = 1; \
    shape_[2] = 1;               strides_[2] = 1; \
  } else if ( array.rank() == 2 ) { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = array.shape(1);  strides_[1] = array.stride(1); \
    shape_[2] = 1;               strides_[2] = 1; \
  } else if ( array.rank() == 3 ) { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = array.shape(1);  strides_[1] = array.stride(1); \
    shape_[2] = array.shape(2);  strides_[2] = array.stride(2); \
  } else { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = array.shape(1);  strides_[1] = array.stride(1); \
    shape_[2] = array.stride(1); strides_[2] = 1; \
  } \
  size_ = shape_[0]*shape_[1]*shape_[2]; \
  ASSERT(size() == array.size()); \
} \
template<> \
ArrayView<DATA_TYPE,4>::ArrayView( const Array& array ) : data_( const_cast<DATA_TYPE*>(array.data<DATA_TYPE>()) ) \
{ \
  ASSERT( array.rank() > 0 ); \
  if( array.rank() == 1  ) { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = 1;               strides_[1] = 1; \
    shape_[2] = 1;               strides_[2] = 1; \
    shape_[3] = 1;               strides_[3] = 1; \
  } else if ( array.rank() == 2 ) { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = array.shape(1);  strides_[1] = array.stride(1); \
    shape_[2] = 1;               strides_[2] = 1; \
    shape_[3] = 1;               strides_[3] = 1; \
  } else if ( array.rank() == 3 ) { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = array.shape(1);  strides_[1] = array.stride(1); \
    shape_[2] = array.shape(2);  strides_[2] = array.stride(2); \
    shape_[3] = 1;               strides_[3] = 1; \
  } else if ( array.rank() == 4 ) { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = array.shape(1);  strides_[1] = array.stride(1); \
    shape_[2] = array.shape(2);  strides_[2] = array.stride(2); \
    shape_[3] = array.shape(3);  strides_[3] = array.stride(3); \
  } else { \
    shape_[0] = array.shape(0);  strides_[0] = array.stride(0); \
    shape_[1] = array.shape(1);  strides_[1] = array.stride(1); \
    shape_[2] = array.shape(2);  strides_[2] = array.stride(2); \
    shape_[3] = array.stride(2); strides_[3] = 1; \
  } \
  size_ = shape_[0]*shape_[1]*shape_[2]*shape_[3];\
  ASSERT(size() == array.size()); \
}\

TEMPLATE_SPECIALIZATION(int);
TEMPLATE_SPECIALIZATION(long);
TEMPLATE_SPECIALIZATION(float);
TEMPLATE_SPECIALIZATION(double);

#undef TEMPLATE_SPECIALIZATION

//------------------------------------------------------------------------------------------------------

} // namespace atlas
