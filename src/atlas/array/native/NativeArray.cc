#include <iostream>
#include "atlas/array/Array.h"
#include "atlas/array/ArrayUtil.h"
#include "atlas/array/MakeView.h"
#include "atlas/array/native/NativeDataStore.h"

namespace atlas {
namespace array {

template<typename Value> Array* Array::create( size_t dim0 ) {
  return new ArrayT<Value>( dim0 );
}
template<typename Value> Array* Array::create( size_t dim0, size_t dim1 ) {
  return new ArrayT<Value>( dim0, dim1 );
}
template<typename Value> Array* Array::create( size_t dim0, size_t dim1, size_t dim2 ) {
  return new ArrayT<Value>( dim0, dim1, dim2 );
}
template<typename Value> Array* Array::create( size_t dim0, size_t dim1, size_t dim2, size_t dim3 ) {
  return new ArrayT<Value>( dim0, dim1, dim2, dim3 );
}
template<typename Value> Array* Array::create( size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4 ) {
  return new ArrayT<Value>( dim0, dim1, dim2, dim3, dim4 );
}
template <typename Value> Array* Array::create( const ArrayShape& shape ) {
  return new ArrayT<Value>( shape );
}
template <typename Value> Array* Array::create( const ArrayShape& shape, const ArrayLayout& layout ) {
  return new ArrayT<Value>( shape,layout );
}
template <typename Value> Array* Array::wrap( Value* data, const ArrayShape& shape ) {
  return new ArrayT<Value>( new native::WrappedDataStore<Value>(data), shape );
}
template <typename Value> Array* Array::wrap( Value* data, const ArraySpec& spec ) {
  return new ArrayT<Value>( new native::WrappedDataStore<Value>(data), spec );
}

Array* Array::create( DataType datatype, const ArrayShape& shape )
{
  switch( datatype.kind() )
  {
    case DataType::KIND_REAL64: return new ArrayT<double>(shape);
    case DataType::KIND_REAL32: return new ArrayT<float>(shape);
    case DataType::KIND_INT32:  return new ArrayT<int>(shape);
    case DataType::KIND_INT64:  return new ArrayT<long>(shape);
    case DataType::KIND_UINT64: return new ArrayT<unsigned long>(shape);
    default:
    {
      std::stringstream err; err << "data kind " << datatype.kind() << " not recognised.";
      throw eckit::BadParameter(err.str(),Here());
    }
  }
  return 0;
}

template <typename Value> ArrayT<Value>::ArrayT(ArrayDataStore* ds, const ArraySpec& spec) {
  data_store_ = std::unique_ptr<ArrayDataStore>(ds);
  spec_ = spec;
}

template <typename Value> ArrayT<Value>::ArrayT(size_t dim0) {
  spec_ = ArraySpec(make_shape(dim0));
  data_store_ = std::unique_ptr<ArrayDataStore>( new native::DataStore<Value>( spec_.size() ) );
}
template <typename Value> ArrayT<Value>::ArrayT(size_t dim0, size_t dim1) {
  spec_ = ArraySpec(make_shape(dim0,dim1));
  data_store_ = std::unique_ptr<ArrayDataStore>( new native::DataStore<Value>( spec_.size() ) );
}
template <typename Value> ArrayT<Value>::ArrayT(size_t dim0, size_t dim1, size_t dim2) {
  spec_ = ArraySpec(make_shape(dim0,dim1,dim2));
  data_store_ = std::unique_ptr<ArrayDataStore>( new native::DataStore<Value>( spec_.size() ) );
}
template <typename Value> ArrayT<Value>::ArrayT(size_t dim0, size_t dim1, size_t dim2, size_t dim3) {
  spec_ = ArraySpec(make_shape(dim0,dim1,dim2,dim3));
  data_store_ = std::unique_ptr<ArrayDataStore>( new native::DataStore<Value>( spec_.size() ) );
}
template <typename Value> ArrayT<Value>::ArrayT(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4) {
  spec_ = ArraySpec(make_shape(dim0,dim1,dim2,dim3,dim4));
  data_store_ = std::unique_ptr<ArrayDataStore>( new native::DataStore<Value>( spec_.size() ) );
}

template <typename Value> ArrayT<Value>::ArrayT(const ArrayShape& shape) {
  ASSERT(shape.size()>0);
  size_t size = 1;
  for( size_t j=0; j<shape.size(); ++j ) size *= shape[j];
  data_store_ = std::unique_ptr<ArrayDataStore>( new native::DataStore<Value>( size ) );
  spec_ = ArraySpec(shape);
}

template <typename Value> ArrayT<Value>::ArrayT(const ArrayShape& shape, const ArrayLayout& layout) {
  spec_ = ArraySpec(shape);
  data_store_ = std::unique_ptr<ArrayDataStore>( new native::DataStore<Value>( spec_.size() ) );
  for( size_t j=0; j<layout.size(); ++j )
    ASSERT( spec_.layout()[j] == layout[j] );
}

template <typename Value> ArrayT<Value>::ArrayT(const ArraySpec& spec) {
  if( not spec.contiguous() )     NOTIMP;
  spec_ = spec;
  data_store_ = std::unique_ptr<ArrayDataStore>( new native::DataStore<Value>( spec_.size() ) );
}

template< typename Value >
void ArrayT<Value>::resize( const ArrayShape& _shape )
{
  Array* resized = Array::create<Value>(_shape);

  Value *resized_data = make_storageview<Value>(*resized).data();
  Value *this_data    = make_storageview<Value>(*this).data();

  for( size_t j=0; j<this->size(); ++j ) {
    resized_data[j] = this_data[j];
  }

  data_store_.swap(resized->data_store_);
  spec_ = resized->spec();

  delete resized;
}

template< typename Value >
void ArrayT<Value>::insert(size_t idx1, size_t size1)
{

  ArrayShape nshape = shape();
  if(idx1 > nshape[0]) {
      throw eckit::BadParameter("can not insert into an array at a position beyond its size", Here());
  }
  nshape[0] += size1;

  Array* resized = Array::create<Value>(nshape);

  Value *resized_data = make_storageview<Value>(*resized).data();
  Value *this_data    = make_storageview<Value>(*this).data();

  size_t insert_begin = idx1*this->stride(0);
  size_t insert_end = insert_begin + size1*this->stride(0);

  size_t c(0);
  for( size_t j=0; j<insert_begin; ++j ) {
    resized_data[j] = this_data[c++];
  }
  for( size_t j=insert_end; j<resized->size(); ++j ) {
    resized_data[j] = this_data[c++];
  }
  assert( c == this->size() );

  data_store_.swap(resized->data_store_);
  spec_ = resized->spec();

  delete resized;
}

template< typename Value >
void ArrayT<Value>::resize(size_t size1) { resize( make_shape(size1) ); }

template< typename Value >
void ArrayT<Value>::resize(size_t size1, size_t size2) { resize( make_shape(size1,size2) ); }

template< typename Value >
void ArrayT<Value>::resize(size_t size1, size_t size2, size_t size3) { resize( make_shape(size1,size2,size3) ); }

template< typename Value >
void ArrayT<Value>::resize(size_t size1, size_t size2, size_t size3, size_t size4) { resize( make_shape(size1,size2,size3,size4) ); }

template< typename Value >
void ArrayT<Value>::resize(size_t size1, size_t size2, size_t size3, size_t size4, size_t size5) { resize( make_shape(size1,size2,size3,size4,size5) ); }

template< typename Value >
void ArrayT<Value>::dump(std::ostream& os) const {
  if( not contiguous() ) NOTIMP;

  Value *data = make_storageview<Value>(*this).data();

  for(size_t i=0; i<size(); ++i) {
      os << data[i] << " ";
      if( (i+1)%10 == 0 )
          os << std::endl;
  }
  os << std::endl;
}

//------------------------------------------------------------------------------

template Array* Array::create<int>(size_t);
template Array* Array::create<long>(size_t);
template Array* Array::create<float>(size_t);
template Array* Array::create<double>(size_t);
template Array* Array::create<long unsigned>(size_t);

template Array* Array::create<int>(size_t,size_t);
template Array* Array::create<long>(size_t,size_t);
template Array* Array::create<float>(size_t,size_t);
template Array* Array::create<double>(size_t,size_t);
template Array* Array::create<long unsigned>(size_t,size_t);

template Array* Array::create<int>(size_t,size_t,size_t);
template Array* Array::create<long>(size_t,size_t,size_t);
template Array* Array::create<float>(size_t,size_t,size_t);
template Array* Array::create<double>(size_t,size_t,size_t);
template Array* Array::create<long unsigned>(size_t,size_t,size_t);

template Array* Array::create<int>(size_t,size_t,size_t,size_t);
template Array* Array::create<long>(size_t,size_t,size_t,size_t);
template Array* Array::create<float>(size_t,size_t,size_t,size_t);
template Array* Array::create<double>(size_t,size_t,size_t,size_t);
template Array* Array::create<long unsigned>(size_t,size_t,size_t,size_t);

template Array* Array::create<int>(size_t,size_t,size_t,size_t,size_t);
template Array* Array::create<long>(size_t,size_t,size_t,size_t,size_t);
template Array* Array::create<float>(size_t,size_t,size_t,size_t,size_t);
template Array* Array::create<double>(size_t,size_t,size_t,size_t,size_t);
template Array* Array::create<long unsigned>(size_t,size_t,size_t,size_t,size_t);

template Array* Array::create<int>(const ArrayShape&);
template Array* Array::create<long>(const ArrayShape&);
template Array* Array::create<float>(const ArrayShape&);
template Array* Array::create<double>(const ArrayShape&);
template Array* Array::create<long unsigned>(const ArrayShape&);

template Array* Array::create<int>(const ArrayShape&, const ArrayLayout&);
template Array* Array::create<long>(const ArrayShape&, const ArrayLayout&);
template Array* Array::create<float>(const ArrayShape&, const ArrayLayout&);
template Array* Array::create<double>(const ArrayShape&, const ArrayLayout&);
template Array* Array::create<long unsigned>(const ArrayShape&, const ArrayLayout&);

template Array* Array::wrap<int>(int*, const ArrayShape&);
template Array* Array::wrap<long>(long*, const ArrayShape&);
template Array* Array::wrap<float>(float*, const ArrayShape&);
template Array* Array::wrap<double>(double*, const ArrayShape&);
template Array* Array::wrap<long unsigned>(long unsigned*, const ArrayShape&);

template Array* Array::wrap<int>(int*, const ArraySpec&);
template Array* Array::wrap<long>(long*, const ArraySpec&);
template Array* Array::wrap<float>(float*, const ArraySpec&);
template Array* Array::wrap<double>(double*, const ArraySpec&);
template Array* Array::wrap<long unsigned>(long unsigned*, const ArraySpec&);

template class ArrayT<int>;
template class ArrayT<long>;
template class ArrayT<float>;
template class ArrayT<double>;
template class ArrayT<unsigned long>;

} // namespace array
} // namespace atlas
