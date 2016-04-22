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
#include "atlas/field/Field.h"
#include "atlas/array/Array.h"
#include "atlas/array/IndexView.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

#define TEMPLATE_SPECIALIZATION( DATA_TYPE ) \
template<>\
IndexView <DATA_TYPE, 1 >::IndexView( const ArrayT<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);  shape_[0]=array.shape(0);\
} \
template<> \
IndexView <DATA_TYPE, 1 >::IndexView( const field::Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);  shape_[0]=field.shape(0); \
} \
template<> \
IndexView<DATA_TYPE,2>::IndexView( const ArrayT<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);       shape_[0]=array.shape(0); \
  strides_[1]=array.stride(1);       shape_[1]=array.shape(1); \
  size_ = shape_[0]*shape_[1]; \
} \
template<> \
IndexView<DATA_TYPE,2>::IndexView( const field::Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);       shape_[0]=field.shape(0); \
  strides_[1]=field.stride(1);       shape_[1]=field.shape(1); \
  size_ = shape_[0]*shape_[1]; \
} \
template<> \
IndexView<DATA_TYPE,3>::IndexView( const ArrayT<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);       shape_[0]=array.shape(0); \
  strides_[1]=array.stride(1);       shape_[1]=array.shape(1); \
  strides_[2]=array.stride(2);       shape_[2]=array.shape(2); \
  size_ = shape_[0]*shape_[1]*shape_[2]; \
} \
template<> \
IndexView<DATA_TYPE,3>::IndexView( const field::Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);       shape_[0]=field.shape(0); \
  strides_[1]=field.stride(1);       shape_[1]=field.shape(1); \
  strides_[2]=field.stride(2);       shape_[2]=field.shape(2); \
  size_ = shape_[0]*shape_[1]*shape_[2]; \
} \
template<> \
IndexView<DATA_TYPE,4>::IndexView( const ArrayT<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);       shape_[0]=array.shape(0); \
  strides_[1]=array.stride(1);       shape_[1]=array.shape(1); \
  strides_[2]=array.stride(2);       shape_[2]=array.shape(2); \
  strides_[3]=array.stride(3);       shape_[3]=array.shape(3); \
  size_ = shape_[0]*shape_[1]*shape_[2]*shape_[3]; \
} \
template<> \
IndexView<DATA_TYPE,4>::IndexView( const field::Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);       shape_[0]=field.shape(0); \
  strides_[1]=field.stride(1);       shape_[1]=field.shape(1); \
  strides_[2]=field.stride(2);       shape_[2]=field.shape(2); \
  strides_[3]=field.stride(3);       shape_[3]=field.shape(3); \
  size_ = shape_[0]*shape_[1]*shape_[2]*shape_[3]; \
}

TEMPLATE_SPECIALIZATION(int);

#undef TEMPLATE_SPECIALIZATION

//------------------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
