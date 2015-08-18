/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>
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
        << " with datatype " << array->datatype().str() << " to "
        << DataType::str<DATA_TYPE>();
    throw eckit::BadCast(msg.str(),Here());
  }
  return const_cast<DATA_TYPE*>(array->data());
}

template< typename DATA_TYPE >
void dump_array_data( const Array<DATA_TYPE>& array, std::ostream& os )
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

template <> const int*    ArrayBase::data<int   >() const { return get_array_data<int   >(*this); }
template <>       int*    ArrayBase::data<int   >()       { return get_array_data<int   >(*this); }
template <> const long*   ArrayBase::data<long  >() const { return get_array_data<long  >(*this); }
template <>       long*   ArrayBase::data<long  >()       { return get_array_data<long  >(*this); }
template <> const float*  ArrayBase::data<float >() const { return get_array_data<float >(*this); }
template <>       float*  ArrayBase::data<float >()       { return get_array_data<float >(*this); }
template <> const double* ArrayBase::data<double>() const { return get_array_data<double>(*this); }
template <>       double* ArrayBase::data<double>()       { return get_array_data<double>(*this); }

template <> void Array<int>::dump(std::ostream& os) const { dump_array_data(*this,os); };
template <> void Array<long>::dump(std::ostream& os) const { dump_array_data(*this,os); };
template <> void Array<float>::dump(std::ostream& os) const { dump_array_data(*this,os); };
template <> void Array<double>::dump(std::ostream& os) const { dump_array_data(*this,os); };


ArrayBase* ArrayBase::create( DataType datatype, const ArrayShape& shape )
{
  switch( datatype.kind() )
  {
    case DataType::KIND_REAL64: return new Array<double>(shape);
    case DataType::KIND_REAL32: return new Array<float>(shape);
    case DataType::KIND_INT32:  return new Array<int>(shape);
    case DataType::KIND_INT64:  return new Array<long>(shape);
    default:
    {
      std::stringstream err; err << "data kind " << datatype.kind() << " not recognised.";
      throw eckit::BadParameter(err.str(),Here());
    }
  }
  return 0;
}


} // namespace atlas
