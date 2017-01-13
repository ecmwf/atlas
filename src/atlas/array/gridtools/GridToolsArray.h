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
#include <iosfwd>
#include <iterator>
#include <type_traits>
#include "eckit/memory/Owned.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/array/ArrayUtil.h"
#include "atlas/array/DataType.h"
#include "atlas/array/ArrayView.h"

#include "atlas/array/gridtools/GridToolsTraits.h"
#include "atlas/array/gridtools/GridToolsDataStoreWrapper.h"
#include "atlas/array/ArrayHelpers.h"
#include "atlas/array/array_fwd.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

template<typename DATA_TYPE>
class ArrayT;

template <unsigned int RANK>
struct array_initializer;

template<unsigned int PartDim>
struct array_initializer_partitioned;


//------------------------------------------------------------------------------

#ifdef GTNS
namespace gridtools {
#endif

class Array : public eckit::Owned {
public:
  static Array* create( array::DataType, const ArrayShape& );
  static Array* create( array::DataType );
  static Array* create( const Array& );

private:

  template<typename Value>
  struct storage_creator {
      template<typename UInt, UInt ... Indices>
      static Array* apply(const ArrayShape& shape, ::gridtools::gt_integer_sequence<UInt, Indices...> ) {
          return Array::create<Value>(shape[Indices]...);
      }

  };

  struct storage_resizer {
      template<typename UInt, UInt ... Indices>
      static void apply(Array& array, const ArrayShape& shape, ::gridtools::gt_integer_sequence<UInt, Indices...> ) {
          return array.resize(shape[Indices]...);
      }
  };

public:

  template <typename Value, typename... UInts,
            typename = ::gridtools::all_integers<UInts...> >
  static Array* create(UInts... dims) {
      return create_with_layout<Value, typename atlas::array::default_layout_t<sizeof...(dims)>::type >(dims...);
  }

  template <typename Value,
            typename Layout,
            typename... UInts,
            typename = ::gridtools::all_integers<UInts...> >
  static Array* create_with_layout(UInts... dims) {
    auto gt_data_store_ptr = create_gt_storage<Value, Layout>(dims...);

    Array* array = new ArrayT<Value>(gt_data_store_ptr);

    array->build_spec(gt_data_store_ptr, dims...);

    return array;
  }


  template<typename Value>
  static Array* create(const ArrayShape& shape)
  {
    return new ArrayT<Value>(shape);
  }

  template <typename Value, template <class> class Storage, typename StorageInfo>
  static Array* wrap_array(::gridtools::data_store< Storage<Value>, StorageInfo> * ds, const ArraySpec& spec) {
    assert(ds);
    Array* array = new ArrayT<Value>(ds);
    array->spec_ = spec;

    return array;
  }

  template <typename T> static Array* wrap(T* data, const ArraySpec& spec) {
    ArrayShape const& shape = spec.shape();
    ArrayStrides const& strides = spec.strides();

    assert(shape.size() > 0);
    switch (shape.size()) {
      case 1:
        return wrap_array( wrap_gt_storage<T, 1>(data, get_array_from_vector<1>(shape), get_array_from_vector<1>(strides)), spec);
      case 2:
        return wrap_array( wrap_gt_storage<T, 2>(data, get_array_from_vector<2>(shape), get_array_from_vector<2>(strides)), spec);
      case 3:
        return wrap_array( wrap_gt_storage<T, 3>(data, get_array_from_vector<3>(shape), get_array_from_vector<3>(strides)), spec);
      case 4:
        return wrap_array( wrap_gt_storage<T, 4>(data, get_array_from_vector<4>(shape), get_array_from_vector<4>(strides)), spec);
      case 5:
        return wrap_array( wrap_gt_storage<T, 5>(data, get_array_from_vector<5>(shape), get_array_from_vector<5>(strides)), spec);
      case 6:
        return wrap_array( wrap_gt_storage<T, 6>(data, get_array_from_vector<6>(shape), get_array_from_vector<6>(strides)), spec);
      case 7:
        return wrap_array( wrap_gt_storage<T, 7>(data, get_array_from_vector<7>(shape), get_array_from_vector<7>(strides)), spec);
      case 8:
        return wrap_array( wrap_gt_storage<T, 8>(data, get_array_from_vector<8>(shape), get_array_from_vector<8>(strides)), spec);
      case 9:
        return wrap_array( wrap_gt_storage<T, 9>(data, get_array_from_vector<9>(shape), get_array_from_vector<9>(strides)), spec);
      default: {
        std::stringstream err;
        err << "shape not recognized";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }

  template <typename T> static Array* wrap(T* data, const ArrayShape& shape) {
    return wrap(data,ArraySpec(shape));
  }

public:

  Array(){}
  Array(const ArraySpec& s) : spec_(s) {}

  template<typename DataStore, typename = typename std::enable_if< ::gridtools::is_data_store<DataStore>::value > >
  Array(DataStore* ds) : data_store_( new DataStoreWrapper<DataStore>(ds)) {}

  std::unique_ptr< DataStoreInterface>& data_store() {
      return data_store_; }

  template<typename DataStore, typename ... Dims>
  void build_spec(DataStore* gt_data_store_ptr, Dims...dims) {
      static_assert((::gridtools::is_data_store<DataStore>::value), "Internal Error: passing a non GT data store");

      auto storage_info_ptr = gt_data_store_ptr->get_storage_info_ptr();
      using Layout = typename DataStore::storage_info_t::Layout;

      using seq =
          ::gridtools::apply_gt_integer_sequence<typename ::gridtools::make_gt_integer_sequence<int, sizeof...(dims)>::type>;

      spec_ = ArraySpec(
        ArrayShape{(unsigned long)dims...},
          seq::template apply<
                        std::vector<unsigned long>,
                        get_stride_component<unsigned long, typename get_pack_size<Dims...>::type>::template get_component>(
                        storage_info_ptr),
          seq::template apply<
                        std::vector<unsigned long>,
                        get_layout_map_component<unsigned long, Layout>::template get_component>()
        );
  }

  virtual void* storage() = 0;
  virtual const void* storage() const = 0;

  virtual array::DataType datatype() const = 0;
  virtual size_t sizeof_data() const = 0;

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

  void insert(size_t idx1, size_t size1);

  template <typename... Coords, typename = ::gridtools::all_integers<Coords...> >
  void resize(Coords... c) {
      if(sizeof...(c) != spec_.rank()){
        std::stringstream err; err << "trying to resize an array of rank " << spec_.rank() << " by dimensions with rank " <<
                                      sizeof...(c) << std::endl;
        throw eckit::BadParameter(err.str(),Here());
      }

      check_dimension_lengths(shape(), c...);

      if(!data_store_->is_on_host()) {
          data_store_->clone_from_device();
      }

      Array* array_resized = create(datatype(), ArrayShape{(unsigned int)c...});

      array_initializer<sizeof...(c)>::apply( *this, *array_resized);
      data_store_.swap(array_resized->data_store());

      spec_ = array_resized->spec();

      //TODO when deleting this if seg fault
      delete array_resized;
  }

  void resize(const ArrayShape& shape)
  {
    assert(shape.size() > 0);
    switch (shape.size()) {
      case 1:
        return storage_resizer::apply(*this, shape, ::gridtools::make_gt_integer_sequence<unsigned int, 1>());
      case 2:
        return storage_resizer::apply(*this, shape, ::gridtools::make_gt_integer_sequence<unsigned int, 2>());
      case 3:
        return storage_resizer::apply(*this, shape, ::gridtools::make_gt_integer_sequence<unsigned int, 3>());
      case 4:
        return storage_resizer::apply(*this, shape, ::gridtools::make_gt_integer_sequence<unsigned int, 4>());
      case 5:
        return storage_resizer::apply(*this, shape, ::gridtools::make_gt_integer_sequence<unsigned int, 5>());
      case 6:
        return storage_resizer::apply(*this, shape, ::gridtools::make_gt_integer_sequence<unsigned int, 6>());
      case 7:
        return storage_resizer::apply(*this, shape, ::gridtools::make_gt_integer_sequence<unsigned int, 7>());
      case 8:
        return storage_resizer::apply(*this, shape, ::gridtools::make_gt_integer_sequence<unsigned int, 8>());
      case 9:
        return storage_resizer::apply(*this, shape, ::gridtools::make_gt_integer_sequence<unsigned int, 9>());
      default: {
        std::stringstream err;
        err << "shape not recognized";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }

  virtual void dump(std::ostream& os) const = 0;

  void clone_to_device() const {
      data_store_->clone_to_device();
  }
  void clone_from_device() const {
      data_store_->clone_from_device();
  }

  bool valid() const {
      return data_store_->valid();
  }
  void sync() const {
      data_store_->sync();
  }
  bool is_on_host() const {
      return data_store_->is_on_host();
  }
  bool is_on_device() const {
      return data_store_->is_on_device();
  }

  void reactivate_device_write_views() const {
      data_store_->reactivate_device_write_views();
  }

  void reactivate_host_write_views() const {
      data_store_->reactivate_host_write_views();
  }

protected:
  std::unique_ptr< DataStoreInterface>  data_store_;

public:
  ArraySpec& spec() {return spec_;}

private: // methods
  ArraySpec spec_;
};

#ifdef GTNS
} // namespace gridtools
#endif

//------------------------------------------------------------------------------

template<typename DATA_TYPE>
class ArrayT : public GTArray  {

public:

  template<typename DataStore, typename = typename std::enable_if<::gridtools::is_data_store<DataStore>::value >::type >
  ArrayT(DataStore* ds): owned_(true), Array(ds) {}

  template <typename... UInts,
            typename = ::gridtools::all_integers<UInts...> >
  ArrayT(UInts... dims) :
    owned_(true)
  {
    create_from_variadic_args(dims...);
  }

  ArrayT(const ArrayShape& shape) :
    owned_(true)
  {
    create_from_shape(shape);
  }

  ArrayT(const ArraySpec& spec) :
    owned_(true)
  {
    if( not spec.default_layout() ) NOTIMP;
    if( not spec.contiguous() )     NOTIMP;
    create_from_shape(spec.shape());
  }

  template <typename... UInts,
            typename = ::gridtools::all_integers<UInts...> >
  void create_from_variadic_args(UInts... dims)
  {
    auto gt_storage = create_gt_storage<DATA_TYPE,typename atlas::array::default_layout_t<sizeof...(dims)>::type>(dims...);
    using data_store_t = typename std::remove_pointer<decltype(gt_storage)>::type;
    build_spec(gt_storage,dims...);
    data_store_ = std::unique_ptr< DataStoreInterface>(new DataStoreWrapper<data_store_t>(gt_storage));
  }

  void create_from_shape(const ArrayShape& shape)
  {
    assert(shape.size() > 0);
    switch (shape.size()) {
      case 1: {
        create_from_variadic_args(shape[0]);
        break;
      }
      case 2: {
        create_from_variadic_args(shape[0],shape[1]);
        break;
      }
      case 3: {
        create_from_variadic_args(shape[0],shape[1],shape[2]);
        break;
      }
      case 4: {
        create_from_variadic_args(shape[0],shape[1],shape[2],shape[3]);
        break;
      }
      case 5: {
        create_from_variadic_args(shape[0],shape[1],shape[2],shape[3],shape[4]);
        break;
      }
      default: {
        std::stringstream err;
        err << "shape not recognized";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }

  virtual array::DataType datatype() const { return array::DataType::create<DATA_TYPE>(); }

  size_t sizeof_data() const {return sizeof(DATA_TYPE);}
  virtual void* storage() { return data_store_->void_data_store();}
  virtual const void* storage() const { return data_store_->void_data_store();}

  virtual void dump(std::ostream& os) const {
    os << "\nWARNING: Array dump not implemented\n";
  }

private:
  bool owned_;
};

//------------------------------------------------------------------------------

} // namespace array
} // namespace atlas

#include "atlas/array/gridtools/GridToolsArray_impl.h"
