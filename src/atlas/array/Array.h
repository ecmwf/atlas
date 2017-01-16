/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <vector>
#include "eckit/memory/Owned.h"
#include "atlas/internals/atlas_defines.h"
#include "atlas/array/ArrayUtil.h"
#include "atlas/array/DataType.h"

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE

namespace atlas {
namespace array {

class Array;

template <typename Value>
class ArrayBackend
{
public:
  static Array* create(size_t dim0);
  static Array* create(size_t dim0, size_t dim1);
  static Array* create(size_t dim0, size_t dim1, size_t dim2);
  static Array* create(size_t dim0, size_t dim1, size_t dim2, size_t dim3);
  static Array* create(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4);
  static Array* create( const ArrayShape& );
  static Array* create( const ArrayShape&, const ArrayLayout& );
  static Array* wrap(Value* data, const ArrayShape& );
  static Array* wrap(Value* data, const ArraySpec& );

public:
  static Array* create(ArrayDataStore*, const ArraySpec&);
};
class ArrayBackendResize;
class ArrayBackendInsert;



class Array : public eckit::Owned
{
public:

  static Array* create( array::DataType, const ArrayShape& );

  static Array* create( array::DataType, const ArrayShape&, const ArrayLayout& );

  template<typename Value> static Array* create(size_t dim0) {
    return ArrayBackend<Value>::create(dim0);
  }

  template<typename Value> static Array* create(size_t dim0, size_t dim1) {
    return ArrayBackend<Value>::create(dim0,dim1);
  }

  template<typename Value> static Array* create(size_t dim0, size_t dim1, size_t dim2) {
    return ArrayBackend<Value>::create(dim0,dim1,dim2);
  }

  template<typename Value> static Array* create(size_t dim0, size_t dim1, size_t dim2, size_t dim3) {
    return ArrayBackend<Value>::create(dim0,dim1,dim2,dim3);
  }

  template<typename Value> static Array* create(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4) {
    return ArrayBackend<Value>::create(dim0,dim1,dim2,dim3,dim4);
  }

  template<typename Value> static Array* create(const ArrayShape& shape) {
    return create( array::DataType::create<Value>(), shape );
  }

  template<typename Value> static Array* create(const ArrayShape& shape, const ArrayLayout& layout) {
    return create( array::DataType::create<Value>(), shape, layout );
  }

  template <typename Value> static Array* wrap(Value* data, const ArrayShape& shape) {
    return ArrayBackend<Value>::wrap(data,shape);
  }

  template <typename Value> static Array* wrap(Value* data, const ArraySpec& spec) {
    return ArrayBackend<Value>::wrap(data,spec);
  }

public:

  Array() {}

  Array(const ArraySpec& spec) : spec_(spec) {}

  size_t bytes() const { return sizeof_data() * size();}

  size_t size() const { return spec_.size(); }

  size_t rank() const { return spec_.rank(); }

  size_t stride(size_t i) const { return spec_.strides()[i]; }

  size_t shape(size_t i) const { return spec_.shape()[i]; }

  const ArrayStrides& strides() const { return spec_.strides(); }

  const ArrayShape& shape() const { return spec_.shape(); }

  const std::vector<int>& shapef() const { return spec_.shapef(); }

  const std::vector<int>& stridesf() const { return spec_.stridesf(); }

  bool contiguous() const { return spec_.contiguous(); }

  bool default_layout() const { return spec_.default_layout(); }

  virtual array::DataType datatype() const = 0;

  virtual size_t sizeof_data() const = 0;

  virtual void resize(const ArrayShape& shape) = 0;
  virtual void resize(size_t dim0) = 0;
  virtual void resize(size_t dim0, size_t dim1) = 0;
  virtual void resize(size_t dim0, size_t dim1, size_t dim2) = 0;
  virtual void resize(size_t dim0, size_t dim1, size_t dim2, size_t dim3) = 0;
  virtual void resize(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4) = 0;

  virtual void insert(size_t idx1, size_t size1) = 0;

  virtual void dump(std::ostream& os) const = 0;

  virtual void* storage() { return data_store_->void_data_store();}

  virtual const void* storage() const { return data_store_->void_data_store();}

  void clone_to_device() const { data_store_->clone_to_device(); }

  void clone_from_device() const { data_store_->clone_from_device(); }

  bool valid() const { return data_store_->valid(); }

  void sync() const { data_store_->sync(); }

  bool is_on_host() const { return data_store_->is_on_host(); }

  bool is_on_device() const { return data_store_->is_on_device(); }

  void reactivate_device_write_views() const { data_store_->reactivate_device_write_views(); }

  void reactivate_host_write_views() const { data_store_->reactivate_host_write_views(); }

  ArraySpec& spec() {return spec_;}

public:
  ArraySpec spec_;
  std::unique_ptr<ArrayDataStore> data_store_;


private:

  friend class ArrayBackendResize;
  friend class ArrayBackendInsert;

};

template<typename DATA_TYPE>
class ArrayT : public Array  {

public:

  // These constructors are used to create an ArrayT on the stack.
  ArrayT(size_t dim0);
  ArrayT(size_t dim0, size_t dim1);
  ArrayT(size_t dim0, size_t dim1, size_t dim2);
  ArrayT(size_t dim0, size_t dim1, size_t dim2, size_t dim3);
  ArrayT(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4);

  ArrayT(const ArrayShape&);
  ArrayT(const ArrayShape&, const ArrayLayout&);

  ArrayT(const ArraySpec&);

public:

  virtual void insert(size_t idx1, size_t size1);

  virtual void resize(const ArrayShape&);
  virtual void resize(size_t dim0);
  virtual void resize(size_t dim0, size_t dim1);
  virtual void resize(size_t dim0, size_t dim1, size_t dim2);
  virtual void resize(size_t dim0, size_t dim1, size_t dim2, size_t dim3);
  virtual void resize(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4);

  virtual array::DataType datatype() const { return array::DataType::create<DATA_TYPE>(); }

  virtual size_t sizeof_data() const {return sizeof(DATA_TYPE);}

  virtual void dump(std::ostream& os) const;

private:
  // This constructor is used through the Array::create() or the Array::wrap() methods
  template <typename T>
  friend class ArrayBackend;
  ArrayT(ArrayDataStore*, const ArraySpec&);
};


} // namespace array
} // namespace atlas

#else
#define GT_FUNCTION
//------------------------------------------------------------------------------

namespace atlas {
namespace array {


template<typename DATA_TYPE>
class ArrayT;

class Array : public eckit::Owned {
public:
  static Array* create( array::DataType, const ArrayShape& );
  static Array* create( array::DataType );
  static Array* create( const Array& );

  template <typename T> static Array* create(const ArrayShape& s);
  template <typename T> static Array* create(size_t size);
  template <typename T> static Array* create(size_t size1, size_t size2);
  template <typename T> static Array* create(size_t size1, size_t size2, size_t size3);
  template <typename T> static Array* create(size_t size1, size_t size2, size_t size3, size_t size4);

  template <typename T> static Array* wrap(T data[], const ArraySpec&);
  template <typename T> static Array* wrap(T data[], const ArrayShape&);

public:


  Array(){}
  Array(const ArraySpec& s) : spec_(s) {}

  virtual void* storage() = 0;
  virtual const void* storage() const = 0;

  virtual array::DataType datatype() const = 0;
  virtual size_t sizeof_data() const = 0;

  size_t bytes() const { return sizeof_data() * size();}

  virtual void dump(std::ostream& os) const = 0;

  void resize(const ArrayShape&);

  void resize(size_t size1);

  void resize(size_t size1, size_t size2);

  void resize(size_t size1, size_t size2, size_t size3);

  void resize(size_t size1, size_t size2, size_t size3, size_t size4);

  void insert(size_t idx1, size_t size1);

  size_t size() const { return spec_.size(); }

  size_t rank() const { return spec_.rank(); }

  size_t stride(size_t i) const { return spec_.strides()[i]; }

  size_t shape(size_t i) const { return spec_.shape()[i]; }

  const ArrayStrides& strides() const { return spec_.strides(); }

  const ArrayShape& shape() const { return spec_.shape(); }

  const std::vector<int>& shapef() const { return spec_.shapef(); }

  const std::vector<int>& stridesf() const { return spec_.stridesf(); }

  bool contiguous() const { return spec_.contiguous(); }

  void clone_to_device() const {
      /* ignore */
  }
  void clone_from_device() const {
      /* ignore */
  }

  bool valid() const {
      return true;
  }
  void sync() const {
      /* ignore */;
  }
  bool is_on_host() const {
      return true;
  }
  bool is_on_device() const {
      return false;
  }

  void reactivate_device_write_views() const {
    /* ignore */
  }

  void reactivate_host_write_views() const {
    /* ignore */
  }

private:
  virtual void resize_data( size_t size )=0;
  virtual void insert_data(size_t idx1, size_t size1)=0;

public:
  ArraySpec& spec() {return spec_;}

private: // methods
  ArraySpec spec_;
};

//------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayT : public Array  {
public:

  typedef typename remove_const<DATA_TYPE>::type  value_type;
  typedef typename add_const<DATA_TYPE>::type     const_value_type;

public:

  ArrayT(): owned_(true) {}

  ArrayT(const ArrayShape& shape): owned_(true)                                { resize(shape); }

  ArrayT(size_t size): owned_(true)                                            { resize( make_shape(size) ); }

  ArrayT(size_t size1, size_t size2): owned_(true)                             { resize( make_shape(size1,size2) ); }

  ArrayT(size_t size1, size_t size2, size_t size3): owned_(true)               { resize( make_shape(size1,size2,size3) ); }

  ArrayT(size_t size1, size_t size2, size_t size3, size_t size4): owned_(true) { resize( make_shape(size1,size2,size3,size4) ); }

  ArrayT(DATA_TYPE data[], const ArraySpec& spec):
    Array(spec),
    owned_(false)
  { wrap(data); }

  ArrayT(DATA_TYPE data[], const ArrayShape& shape):
    Array(ArraySpec(shape)),
    owned_(false)
  { wrap(data); }

public:

  virtual array::DataType datatype() const { return array::DataType::create<DATA_TYPE>(); }

  size_t sizeof_data() const {return sizeof(DATA_TYPE);}

  virtual void dump(std::ostream& os) const;

  virtual void* storage()             { return (data_); }
  virtual const void* storage() const { return (data_); }

private:

  virtual void resize_data( size_t size );
  virtual void insert_data(size_t idx1, size_t size1);
  void wrap(DATA_TYPE data[]);

private:
  bool owned_;
  std::vector<DATA_TYPE> owned_data_;
  DATA_TYPE* data_;
};


template< typename DATA_TYPE>
void ArrayT<DATA_TYPE>::resize_data( size_t size )
{
  if( !owned_ ) throw eckit::SeriousBug("Cannot resize data that is not owned");
  owned_data_.resize( size );
  data_ = owned_data_.data();
}

template< typename DATA_TYPE>
void ArrayT<DATA_TYPE>::insert_data(size_t pos, size_t size)
{
  if( !owned_ ) throw eckit::SeriousBug("Cannot resize data that is not owned");
  owned_data_.insert(owned_data_.begin()+pos,size,0);
  data_ = owned_data_.data();
}


template <typename DATA_TYPE>
void ArrayT<DATA_TYPE>::wrap(DATA_TYPE data[])
{
  data_ = data;
}

//------------------------------------------------------------------------------

template <typename T> Array* Array::create(const ArrayShape& s)
{ return create(array::DataType::create<T>(),s); }

template <typename T> Array* Array::create(size_t size)
{ return create(array::DataType::create<T>(),make_shape(size)); }

template <typename T> Array* Array::create(size_t size1, size_t size2)
{ return create(array::DataType::create<T>(),make_shape(size1,size2)); }

template <typename T> Array* Array::create(size_t size1, size_t size2, size_t size3)
{ return create(array::DataType::create<T>(),make_shape(size1,size2,size3)); }

template <typename T> Array* Array::create(size_t size1, size_t size2, size_t size3, size_t size4)
{ return create(array::DataType::create<T>(),make_shape(size1,size2,size3,size4)); }

} // namespace array
} // namespace atlas

#endif
