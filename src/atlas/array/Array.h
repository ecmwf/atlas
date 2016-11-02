/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Array_h
#define atlas_Array_h

#include <vector>
#include <iosfwd>
#include <iterator>
#include "eckit/memory/Owned.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/array/ArrayUtil.h"
#include "atlas/array/DataType.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/GridToolsTraits.h"

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
#include "GridToolsDataStoreWrapper.h"
#endif

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE

template<typename DATA_TYPE>
class ArrayT;

#endif

class Array : public eckit::Owned {
public:
  static Array* create( array::DataType, const ArrayShape& );
  static Array* create( array::DataType );
  static Array* create( const Array& );

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE

private:

  //indirection around C++11 sizeof... since it is buggy for nvcc and cray
  template<typename ...T>
  struct get_pack_size {
      using type = gridtools::static_uint< sizeof...(T) >;
  };

  template<typename Value>
  struct storage_creator {
      template<typename UInt, UInt ... Indices>
      static Array* apply(const ArrayShape& shape, gridtools::gt_integer_sequence<UInt, Indices...> ) {
          return Array::create<Value>(shape[Indices]...);
      }

  };

public:
  template <typename Value, typename... UInts>
  static gridtools::storage_traits<BACKEND>::data_store_t<
      Value, gridtools::storage_traits<BACKEND>::storage_info_t<0, get_pack_size<UInts...>::type::value> >* create_storage_(UInts... dims) {
    static_assert((sizeof...(dims) > 0), "Error: can not create storages without any dimension");

    typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, get_pack_size<UInts...>::type::value> storage_info_ty;
    storage_info_ty si(dims...);

    typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;
    data_store_t* ds = new data_store_t(si);
    ds->allocate();

    return ds;
  }

public:

  template<typename Value, typename ... UInts, typename = gridtools::all_integers<UInts...> >
  static Array* create(UInts... dims)
  {
      Array* array = new ArrayT<Value>(create_storage_<Value>(dims...));
      array->set_spec(dims...);
      return array;
  }

  template<typename Value>
  static Array* create(const ArrayShape& shape)
  {
    assert(shape.size() > 0);
    switch (shape.size()) {
      case 1:
        return storage_creator<Value>::apply(shape, gridtools::make_gt_integer_sequence<unsigned int, 1>());
      case 2:
        return storage_creator<Value>::apply(shape, gridtools::make_gt_integer_sequence<unsigned int, 2>());
      case 3:
        return storage_creator<Value>::apply(shape, gridtools::make_gt_integer_sequence<unsigned int, 3>());
      case 4:
        return storage_creator<Value>::apply(shape, gridtools::make_gt_integer_sequence<unsigned int, 4>());
      case 5:
        return storage_creator<Value>::apply(shape, gridtools::make_gt_integer_sequence<unsigned int, 5>());
      case 6:
        return storage_creator<Value>::apply(shape, gridtools::make_gt_integer_sequence<unsigned int, 6>());
      case 7:
        return storage_creator<Value>::apply(shape, gridtools::make_gt_integer_sequence<unsigned int, 7>());
      case 8:
        return storage_creator<Value>::apply(shape, gridtools::make_gt_integer_sequence<unsigned int, 8>());
      case 9:
        return storage_creator<Value>::apply(shape, gridtools::make_gt_integer_sequence<unsigned int, 9>());
      default: {
        std::stringstream err;
        err << "shape not recognized";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }

#else

  template <typename T> static Array* create(const ArrayShape& s);
  template <typename T> static Array* create(size_t size);
  template <typename T> static Array* create(size_t size1, size_t size2);
  template <typename T> static Array* create(size_t size1, size_t size2, size_t size3);
  template <typename T> static Array* create(size_t size1, size_t size2, size_t size3, size_t size4);

#endif

  template <typename T> static Array* create();
  template <typename T> static Array* wrap(T data[], const ArraySpec&);
  template <typename T> static Array* wrap(T data[], const ArrayShape&);

public:

  Array(){}
  Array(const ArraySpec& s) : spec_(s) {}
#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
  template<typename DataStore, typename = typename std::enable_if< gridtools::is_data_store<DataStore>::value > >
  Array(DataStore* ds) : data_store_( new DataStoreWrapper<DataStore>(ds)) {}

#endif

  virtual array::DataType datatype() const = 0;
#ifndef ATLAS_HAVE_GRIDTOOLS_STORAGE
  virtual double bytes() const = 0;
  virtual void dump(std::ostream& os) const = 0;

  void resize(const ArrayShape&);

  void resize(size_t size1);

  void resize(size_t size1, size_t size2);

  void resize(size_t size1, size_t size2, size_t size3);

  void resize(size_t size1, size_t size2, size_t size3, size_t size4);

  void insert(size_t idx1, size_t size1);
#endif

  size_t size() const { return spec_.size(); }

  size_t rank() const { return spec_.rank(); }

  size_t stride(size_t i) const { return spec_.strides()[i]; }

  size_t shape(size_t i) const { return spec_.shape()[i]; }

  const ArrayStrides& strides() const { return spec_.strides(); }

  const ArrayShape& shape() const { return spec_.shape(); }

  const std::vector<int>& shapef() const { return spec_.shapef(); }

  const std::vector<int>& stridesf() const { return spec_.stridesf(); }

  bool contiguous() const { return spec_.contiguous(); }

#ifndef ATLAS_HAVE_GRIDTOOLS_STORAGE
  void operator=( const Array &array ) { return assign(array); }

  virtual void assign( const Array& )=0;

  virtual void resize_data( size_t size )=0;
  virtual void insert_data(size_t idx1, size_t size1)=0;

  void operator=( const Array &array ) { return assign(array); }
#endif

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
  template<typename ... UInt>
  void set_spec(UInt... dims)
  {
      spec_ = ArraySpec(ArrayShape{dims...});
  }

  void clone_to_device() const {
      data_store_->clone_to_device();
  }
  void clone_from_device() const {
      data_store_->clone_from_device();
  }

  bool valid() const {
      data_store_->valid();
  }
  void sync() const {
      data_store_->sync();
  }
  bool is_on_host() const {
      data_store_->is_on_host();
  }
  bool is_on_device() const {
      data_store_->is_on_device();
  }

  void reactivate_device_write_views() const {
      data_store_->reactivate_device_write_views();
  }

  void reactivate_host_write_views() const {
      data_store_->reactivate_host_write_views();
  }

private:
  std::unique_ptr< DataStoreInterface>  data_store_;

#endif

private: // methods
  ArraySpec spec_;
};

//------------------------------------------------------------------------------

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
template<typename DATA_TYPE>
class ArrayT : public Array  {
public:

public:

  template<typename DataStore, typename = typename std::enable_if<gridtools::is_data_store<DataStore>::value >::type >
  ArrayT(DataStore* ds): owned_(true), data_(static_cast<void*>(ds)), Array(ds) {}

  virtual array::DataType datatype() const { return array::DataType::create<DATA_TYPE>(); }

  void* data() { return data_;}
private:
  bool owned_;
  void* data_;
};

#else
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
  virtual double bytes() const { return sizeof(DATA_TYPE)*size(); }
  virtual void dump(std::ostream& os) const;

  const DATA_TYPE& operator[](size_t i) const { return *(data()+i); }
        DATA_TYPE& operator[](size_t i)       { return *(data()+i); }

  const DATA_TYPE* data() const { return (data_); }
        DATA_TYPE* data()       { return (data_); }

  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<size(); ++n) data()[n]=scalar; }

  virtual void assign( const Array& );

  template< typename RandomAccessIterator >
  void assign( RandomAccessIterator begin, RandomAccessIterator end );

  virtual void insert_data(size_t idx1, size_t size1);

private:

  virtual void resize_data( size_t size );
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

template< typename DATA_TYPE>
template< typename RandomAccessIterator >
void ArrayT<DATA_TYPE>::assign( RandomAccessIterator begin, RandomAccessIterator end )
{
  if( not contiguous() ) NOTIMP;
  if( std::distance(begin,end) != size() ) {
    throw eckit::SeriousBug("Size doesn't match");
  }
  RandomAccessIterator it = begin;
  for( size_t j=0; j<size(); ++j, ++it ) {
    data()[j] = *it;
  }
}

template< typename DATA_TYPE>
void ArrayT<DATA_TYPE>::assign( const Array& other )
{
  if( not contiguous() or not other.contiguous()) NOTIMP;
  resize( other.shape() );
  ASSERT( datatype().kind() == other.datatype().kind() );
  const ArrayT<DATA_TYPE>& other_array = *dynamic_cast< const ArrayT<DATA_TYPE>* >(&other);
  
  const DATA_TYPE* other_data = other_array.data();
  for( size_t j=0; j<size(); ++j )
    data()[j] = other_data[j];
}
#endif

//------------------------------------------------------------------------------


#ifndef ATLAS_HAVE_GRIDTOOLS_STORAGE

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

#endif


} // namespace array
} // namespace atlas

#endif
