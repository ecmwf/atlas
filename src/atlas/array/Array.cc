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
#include "atlas/array/Array_impl.h"

using atlas::array::DataType;

namespace atlas {
namespace array {

#ifndef ATLAS_HAVE_GRIDTOOLS_STORAGE

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
#endif

#ifndef ATLAS_HAVE_GRIDTOOLS_STORAGE

namespace {
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
template <> void ArrayT<int>::dump(std::ostream& os) const { dump_array_data(*this,os); };
template <> void ArrayT<long>::dump(std::ostream& os) const { dump_array_data(*this,os); };
template <> void ArrayT<float>::dump(std::ostream& os) const { dump_array_data(*this,os); };
template <> void ArrayT<double>::dump(std::ostream& os) const { dump_array_data(*this,os); };

#endif

#ifndef ATLAS_HAVE_GRIDTOOLS_STORAGE

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
#else

Array* Array::create( DataType datatype, const ArrayShape& shape )
{
  switch( datatype.kind() )
  {
    case DataType::KIND_REAL64: return create<double>(shape);
    case DataType::KIND_REAL32: return create<float>(shape);
    case DataType::KIND_INT32:  return create<int>(shape);
    case DataType::KIND_INT64:  return create<long>(shape);
    default:
    {
      std::stringstream err; err << "data kind " << datatype.kind() << " not recognised.";
      throw eckit::BadParameter(err.str(),Here());
    }
  }
  return 0;
}

#endif

#ifndef ATLAS_HAVE_GRIDTOOLS_STORAGE

Array* Array::create( const Array& other )
{
  Array* array = Array::create(other.datatype());
  array->assign(other);
  return array;
}

template <> Array* Array::wrap(int data[], const ArraySpec& s) { return new ArrayT<int>(data,s); }
template <> Array* Array::wrap(long data[], const ArraySpec& s) { return new ArrayT<long>(data,s); }
template <> Array* Array::wrap(float data[], const ArraySpec& s) { return new ArrayT<float>(data,s); }
template <> Array* Array::wrap(double data[], const ArraySpec& s) { return new ArrayT<double>(data,s); }

template <> Array* Array::wrap(int data[], const ArrayShape& s) { return new ArrayT<int>(data,s); }
template <> Array* Array::wrap(long data[], const ArrayShape& s) { return new ArrayT<long>(data,s); }
template <> Array* Array::wrap(float data[], const ArrayShape& s) { return new ArrayT<float>(data,s); }
template <> Array* Array::wrap(double data[], const ArrayShape& s) { return new ArrayT<double>(data,s); }
#endif

} // namespace array
} // namespace atlas
