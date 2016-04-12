/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>
#include "atlas/array/Array.h"

using atlas::array::DataType;

namespace atlas {
namespace array {

void Array::resize(const ArrayShape& _shape)
{
  spec_ = ArraySpec(_shape);
  resize_data(spec_.size());
}

void Array::insert(size_t idx1, size_t size1)
{
  size_t old_size = shape(0);
  ArrayShape _shape = shape();
  _shape[0] += size1;
  spec_ = ArraySpec(_shape);
  if( idx1 == old_size ) {
    resize_data(spec_.size());
  }
  else {
    insert_data(idx1*stride(0),size1*stride(0));
  }
}


void Array::resize(size_t size1) { resize( make_shape(size1) ); }

void Array::resize(size_t size1, size_t size2) { resize( make_shape(size1,size2) ); }

void Array::resize(size_t size1, size_t size2, size_t size3) { resize( make_shape(size1,size2,size3) ); }

void Array::resize(size_t size1, size_t size2, size_t size3, size_t size4) { resize( make_shape(size1,size2,size3,size4) ); }

namespace {
template< typename DATA_TYPE >
DATA_TYPE* get_array_data( const Array& arraybase )
{
  const ArrayT<DATA_TYPE>* array = dynamic_cast< const ArrayT<DATA_TYPE>* >(&arraybase);
  if( array == NULL )
  {
    std::stringstream msg;
    msg << "Could not cast Array "
        << " with datatype " << arraybase.datatype().str() << " to "
        << DataType::str<DATA_TYPE>();
    throw eckit::BadCast(msg.str(),Here());
  }
  return const_cast<DATA_TYPE*>(array->data());
}

template< typename DATA_TYPE >
void dump_array_data( const ArrayT<DATA_TYPE>& array, std::ostream& os )
{
  const DATA_TYPE* data = array.data();
  for(size_t i = 0; i < array.size(); ++i)
  {
      os << data[i] << " ";
      if( (i+1)%10 == 0 )
          os << std::endl;
  }
  os << std::endl;
}

}

template <> const int*    Array::data<int   >() const { return get_array_data<int   >(*this); }
template <>       int*    Array::data<int   >()       { return get_array_data<int   >(*this); }
template <> const long*   Array::data<long  >() const { return get_array_data<long  >(*this); }
template <>       long*   Array::data<long  >()       { return get_array_data<long  >(*this); }
template <> const float*  Array::data<float >() const { return get_array_data<float >(*this); }
template <>       float*  Array::data<float >()       { return get_array_data<float >(*this); }
template <> const double* Array::data<double>() const { return get_array_data<double>(*this); }
template <>       double* Array::data<double>()       { return get_array_data<double>(*this); }

template <> void ArrayT<int>::dump(std::ostream& os) const { dump_array_data(*this,os); };
template <> void ArrayT<long>::dump(std::ostream& os) const { dump_array_data(*this,os); };
template <> void ArrayT<float>::dump(std::ostream& os) const { dump_array_data(*this,os); };
template <> void ArrayT<double>::dump(std::ostream& os) const { dump_array_data(*this,os); };


Array* Array::create( DataType datatype, const ArrayShape& shape )
{
  switch( datatype.kind() )
  {
    case DataType::KIND_REAL64: return new ArrayT<double>(shape);
    case DataType::KIND_REAL32: return new ArrayT<float>(shape);
    case DataType::KIND_INT32:  return new ArrayT<int>(shape);
    case DataType::KIND_INT64:  return new ArrayT<long>(shape);
    default:
    {
      std::stringstream err; err << "data kind " << datatype.kind() << " not recognised.";
      throw eckit::BadParameter(err.str(),Here());
    }
  }
  return 0;
}

Array* Array::create( DataType datatype )
{
  switch( datatype.kind() )
  {
    case DataType::KIND_REAL64: return new ArrayT<double>();
    case DataType::KIND_REAL32: return new ArrayT<float>();
    case DataType::KIND_INT32:  return new ArrayT<int>();
    case DataType::KIND_INT64:  return new ArrayT<long>();
    default:
    {
      std::stringstream err; err << "data kind " << datatype.kind() << " not recognised.";
      throw eckit::BadParameter(err.str(),Here());
    }
  }
  return 0;
}

Array* Array::create( const Array& other )
{
  Array* array = Array::create(other.datatype());
  array->assign(other);
  return array;
}

template <> Array* Array::wrap(int data[], const ArrayShape& s) { return new ArrayT<int>(data,s); }
template <> Array* Array::wrap(long data[], const ArrayShape& s) { return new ArrayT<long>(data,s); }
template <> Array* Array::wrap(float data[], const ArrayShape& s) { return new ArrayT<float>(data,s); }
template <> Array* Array::wrap(double data[], const ArrayShape& s) { return new ArrayT<double>(data,s); }

} // namespace array
} // namespace atlas
