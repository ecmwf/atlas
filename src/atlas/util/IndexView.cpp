/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <stdexcept>

#include "atlas/mesh/Field.hpp"
#include "atlas/util/Array.hpp"
#include "atlas/util/IndexView.hpp"

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

#define TEMPLATE_SPECIALIZATION( DATA_TYPE ) \
template<>\
IndexView <DATA_TYPE, 1 >::IndexView( const Array<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);  extents_[0]=array.extent(0);\
} \
template<> \
IndexView <DATA_TYPE, 1 >::IndexView( const Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);  extents_[0]=field.extent(0); \
} \
template<> \
IndexView<DATA_TYPE,2>::IndexView( const Array<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);       extents_[0]=array.extent(0); \
  strides_[1]=array.stride(1);       extents_[1]=array.extent(1); \
} \
template<> \
IndexView<DATA_TYPE,2>::IndexView( const Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);       extents_[0]=field.extent(0); \
  strides_[1]=field.stride(1);       extents_[1]=field.extent(1); \
} \
template<> \
IndexView<DATA_TYPE,3>::IndexView( const Array<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);       extents_[0]=array.extent(0); \
  strides_[1]=array.stride(1);       extents_[1]=array.extent(1); \
  strides_[2]=array.stride(2);       extents_[2]=array.extent(2); \
} \
template<> \
IndexView<DATA_TYPE,3>::IndexView( const Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);       extents_[0]=field.extent(0); \
  strides_[1]=field.stride(1);       extents_[1]=field.extent(1); \
  strides_[2]=field.stride(2);       extents_[2]=field.extent(2); \
} \
template<> \
IndexView<DATA_TYPE,4>::IndexView( const Array<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);       extents_[0]=array.extent(0); \
  strides_[1]=array.stride(1);       extents_[1]=array.extent(1); \
  strides_[2]=array.stride(2);       extents_[2]=array.extent(2); \
  strides_[3]=array.stride(3);       extents_[3]=array.extent(3); \
} \
template<> \
IndexView<DATA_TYPE,4>::IndexView( const Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);       extents_[0]=field.extent(0); \
  strides_[1]=field.stride(1);       extents_[1]=field.extent(1); \
  strides_[2]=field.stride(2);       extents_[2]=field.extent(2); \
  strides_[3]=field.stride(3);       extents_[3]=field.extent(3); \
}

TEMPLATE_SPECIALIZATION(int);

#undef TEMPLATE_SPECIALIZATION

//------------------------------------------------------------------------------------------------------

} // namespace atlas
