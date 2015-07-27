/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <vector>
#include "atlas/util/Array.h"

namespace atlas {

void ArrayBase::resize(const ArrayShape& _shape)
{
  spec_ = ArraySpec(_shape);
  resize_data(spec_.size());
}


void ArrayBase::resize(size_t size1) { resize( make_shape(size1) ); }

void ArrayBase::resize(size_t size1, size_t size2) { resize( make_shape(size1,size2) ); }

void ArrayBase::resize(size_t size1, size_t size2, size_t size3) { resize( make_shape(size1,size2,size3) ); }

void ArrayBase::resize(size_t size1, size_t size2, size_t size3, size_t size4) { resize( make_shape(size1,size2,size3,size4) ); }

namespace {
template< typename DATA_TYPE >
DATA_TYPE* get_array_data( const ArrayBase& arraybase )
{
  const Array<DATA_TYPE>* array = dynamic_cast< const Array<DATA_TYPE>* >(&arraybase);
  if( array == NULL )
  {
    std::stringstream msg;
    msg << "Could not cast Array "
        << " with datatype " << array->datatype() << " to "
        << DataType::datatype<DATA_TYPE>();
    throw eckit::BadCast(msg.str(),Here());
  }
  return const_cast<DATA_TYPE*>(array->data());
}
}

template <> const int*    ArrayBase::data<int   >() const { return get_array_data<int   >(*this); }
template <>       int*    ArrayBase::data<int   >()       { return get_array_data<int   >(*this); }
template <> const long*   ArrayBase::data<long  >() const { return get_array_data<long  >(*this); }
template <>       long*   ArrayBase::data<long  >()       { return get_array_data<long  >(*this); }
template <> const float*  ArrayBase::data<float >() const { return get_array_data<float >(*this); }
template <>       float*  ArrayBase::data<float >()       { return get_array_data<float >(*this); }
template <> const double* ArrayBase::data<double>() const { return get_array_data<double>(*this); }
template <>       double* ArrayBase::data<double>()       { return get_array_data<double>(*this); }


} // namespace atlas
