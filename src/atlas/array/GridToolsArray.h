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

#include "atlas/array/GridToolsTraits.h"
#include "atlas/array/GridToolsDataStoreWrapper.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

template<typename DATA_TYPE>
class ArrayT;

template <unsigned int RANK>
struct array_initializer;

template<unsigned int PartDim>
struct array_initializer_partitioned;

template<unsigned int NDims>
struct default_layout {

    template<typename T>
    struct get_layout;

    template<typename UInt, UInt ... Indices>
    struct get_layout<gridtools::gt_integer_sequence<UInt, Indices...> >
    {
        using type = gridtools::layout_map<Indices...>;
    };

    using type = typename get_layout< typename gridtools::make_gt_integer_sequence<unsigned int, NDims>::type >::type;
};

template <unsigned int NDims>
std::array<unsigned int, NDims> get_array_from_vector(std::vector<size_t> const& values) {
    std::array<unsigned int, NDims> array;
    std::copy(values.begin(), values.end(), array.begin());
    return array;
}

template <unsigned int TotalDims, unsigned int Dim, typename = void>
struct check_dimension_lengths_impl {
  template <typename FirstDim, typename... Dims>
  static void apply(ArrayShape const& shape, FirstDim first_dim, Dims... d) {
    if (first_dim < shape[Dim]) {
      std::stringstream err;
      err << "Attempt to resize array with original size for dimension " << Dim << " of " << shape[Dim] << " by "
          << first_dim << std::endl;
      throw eckit::BadParameter(err.str(), Here());
    }
    check_dimension_lengths_impl<TotalDims, Dim + 1>::apply(shape, d...);
  }
};

template <unsigned int TotalDims, unsigned int Dim>
struct check_dimension_lengths_impl<TotalDims, Dim, typename std::enable_if< (Dim == TotalDims-1)>::type > {
  template <typename FirstDim>
  static void apply(ArrayShape const& shape, FirstDim first_dim) {
    if (first_dim < shape[Dim]) {
      std::stringstream err;
      err << "Attempt to resize array with original size for dimension " << Dim - 1 << " of "
          << shape[Dim - 1] << " by " << first_dim << std::endl;
      throw eckit::BadParameter(err.str(), Here());
    }
  }
};

template<typename ... Dims>
void check_dimension_lengths(ArrayShape const&  shape, Dims...d) {
    check_dimension_lengths_impl<sizeof...(d), 0>::apply(shape, d...);
}

class Array : public eckit::Owned {
public:
  static Array* create( array::DataType, const ArrayShape& );
  static Array* create( array::DataType );
  static Array* create( const Array& );


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

  struct storage_resizer {
      template<typename UInt, UInt ... Indices>
      static void apply(Array& array, const ArrayShape& shape, gridtools::gt_integer_sequence<UInt, Indices...> ) {
          return array.resize(shape[Indices]...);
      }
  };

public:
  template <typename Value, typename LayoutMap, typename... UInts>
  static gridtools::storage_traits<BACKEND>::data_store_t<
      Value, gridtools::storage_traits<BACKEND>::storage_info_t<
                 0, get_pack_size<UInts...>::type::value,
                 typename gridtools::zero_halo<get_pack_size<UInts...>::type::value>::type, LayoutMap> >*
      create_storage_(UInts... dims) {
    static_assert((sizeof...(dims) > 0), "Error: can not create storages without any dimension");

    constexpr static unsigned int ndims = get_pack_size<UInts...>::type::value;
    typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, ndims, typename gridtools::zero_halo<ndims>::type,
                                                               LayoutMap> storage_info_ty;
    storage_info_ty si(dims...);

    typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;
    data_store_t* ds = new data_store_t(si);
    ds->allocate();

    return ds;
  }

  template <typename Value, unsigned int NDims>
  static gridtools::storage_traits<BACKEND>::data_store_t<
      Value, gridtools::storage_traits<BACKEND>::storage_info_t<
                 0, NDims,
                 typename gridtools::zero_halo<NDims>::type,
                 typename default_layout<NDims>::type > >* wrap_storage_(Value* data,
          std::array<unsigned int, NDims>&& shape, std::array<unsigned int, NDims>&& strides) {

    static_assert((NDims > 0), "Error: can not create storages without any dimension");
    typedef gridtools::storage_traits<BACKEND>::storage_info_t<
        0, NDims, typename gridtools::zero_halo<NDims>::type,
        typename default_layout<NDims>::type> storage_info_ty;
    storage_info_ty si(shape, strides);

    typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;
    data_store_t* ds = new data_store_t(si, data);

    return ds;
  }

public:

  template <typename Value, typename LayoutMap>
  struct get_layout_map_component {
    // TODO: static_assert( gridtools::is_layout_map<LayoutMap>(), "Error: not a layout_map" );
    template <int Idx>
    struct get_component {
      GT_FUNCTION
      constexpr get_component() {}

      GT_FUNCTION constexpr static Value apply() {
        return LayoutMap::template at<Idx>();
      }
    };
  };

  template <typename Value, typename RANKS>
  struct get_stride_component {
    template <int Idx>
    struct get_component {
      GT_FUNCTION
      constexpr get_component() {}

      template <typename StorageInfoPtr>
      GT_FUNCTION constexpr static Value apply(StorageInfoPtr a) {
        static_assert((gridtools::is_storage_info<typename std::remove_pointer<StorageInfoPtr>::type>::value),
                      "Error: not a storage_info");
        return a->template stride<Idx>();
      }
    };
  };

  template < int Idx >
  struct get_shape_component {

      GT_FUNCTION
      constexpr get_shape_component() {}

      template < typename StorageInfoPtr>
      GT_FUNCTION constexpr static int apply(StorageInfoPtr a) {
          static_assert((gridtools::is_storage_info<typename std::remove_pointer<StorageInfoPtr>::type >::value ), "Error: not a storage_info");
          return a->template dim<Idx>();
      }
  };

  template <typename Layout>
  struct get_shapef_component {
    template <int Idx>
    struct get_component {
      GT_FUNCTION
      constexpr get_component() {}

      template <typename StorageInfoPtr>
      GT_FUNCTION constexpr static int apply(StorageInfoPtr a) {
        static_assert((gridtools::is_storage_info<typename std::remove_pointer<StorageInfoPtr>::type>::value),
                      "Error: not a storage_info");
        return a->template dim<Layout::template at<Layout::length - Idx - 1>() >();
      }
    };
  };

  template <typename Value, typename... UInts,
            typename = gridtools::all_integers<UInts...> >
  static Array* create(UInts... dims) {
      return create_with_layout<Value, typename default_layout<sizeof...(dims)>::type >(dims...);
  }

  template <typename Value,
            typename Layout,
            typename... UInts,
            typename = gridtools::all_integers<UInts...> >
  static Array* create_with_layout(UInts... dims) {
    auto gt_data_store_ptr = create_storage_<Value, Layout>(dims...);

    Array* array = new ArrayT<Value>(gt_data_store_ptr);

    array->build_spec(gt_data_store_ptr, dims...);

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

  template <typename T> static Array* create();

  template <typename Value, template <class> class Storage, typename StorageInfo>
  static Array* wrap_array(gridtools::data_store< Storage<Value>, StorageInfo> * ds, const ArraySpec& spec) {
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
        return wrap_array( wrap_storage_<T, 1>(data, get_array_from_vector<1>(strides), get_array_from_vector<1>(strides)), spec);
      case 2:
        return wrap_array( wrap_storage_<T, 2>(data, get_array_from_vector<2>(strides), get_array_from_vector<2>(strides)), spec);
      case 3:
        return wrap_array( wrap_storage_<T, 3>(data, get_array_from_vector<3>(strides), get_array_from_vector<3>(strides)), spec);
      case 4:
        return wrap_array( wrap_storage_<T, 4>(data, get_array_from_vector<4>(strides), get_array_from_vector<4>(strides)), spec);
      case 5:
        return wrap_array( wrap_storage_<T, 5>(data, get_array_from_vector<5>(strides), get_array_from_vector<5>(strides)), spec);
      case 6:
        return wrap_array( wrap_storage_<T, 6>(data, get_array_from_vector<6>(strides), get_array_from_vector<6>(strides)), spec);
      case 7:
        return wrap_array( wrap_storage_<T, 7>(data, get_array_from_vector<7>(strides), get_array_from_vector<7>(strides)), spec);
      case 8:
        return wrap_array( wrap_storage_<T, 8>(data, get_array_from_vector<8>(strides), get_array_from_vector<8>(strides)), spec);
      case 9:
        return wrap_array( wrap_storage_<T, 9>(data, get_array_from_vector<9>(strides), get_array_from_vector<9>(strides)), spec);
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

  template<typename DataStore, typename = typename std::enable_if< gridtools::is_data_store<DataStore>::value > >
  Array(DataStore* ds) : data_store_( new DataStoreWrapper<DataStore>(ds)) {}

  std::unique_ptr< DataStoreInterface>& data_store() { return data_store_; }

  template<typename DataStore, typename ... Dims>
  void build_spec(DataStore* gt_data_store_ptr, Dims...dims) {
      static_assert((gridtools::is_data_store<DataStore>::value), "Internal Error: passing a non GT data store");

      auto storage_info_ptr = gt_data_store_ptr->get_storage_info_ptr();
      using Layout = typename DataStore::storage_info_t::Layout;

      using seq =
          gridtools::apply_gt_integer_sequence<typename gridtools::make_gt_integer_sequence<int, sizeof...(dims)>::type>;

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

  template <typename... Coords, typename = gridtools::all_integers<Coords...> >
  GT_FUNCTION
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

      Array* array_resized = Array::create(datatype(), ArrayShape{(unsigned int)c...});

      array_initializer<sizeof...(c)>::apply( *this, *array_resized);
      data_store_.swap(array_resized->data_store());

      spec_ = array_resized->spec();

      //TODO when deleting this if seg fault
//      delete array_resized;
  }

  void resize(const ArrayShape& shape)
  {
    assert(shape.size() > 0);
    switch (shape.size()) {
      case 1:
        return storage_resizer::apply(*this, shape, gridtools::make_gt_integer_sequence<unsigned int, 1>());
      case 2:
        return storage_resizer::apply(*this, shape, gridtools::make_gt_integer_sequence<unsigned int, 2>());
      case 3:
        return storage_resizer::apply(*this, shape, gridtools::make_gt_integer_sequence<unsigned int, 3>());
      case 4:
        return storage_resizer::apply(*this, shape, gridtools::make_gt_integer_sequence<unsigned int, 4>());
      case 5:
        return storage_resizer::apply(*this, shape, gridtools::make_gt_integer_sequence<unsigned int, 5>());
      case 6:
        return storage_resizer::apply(*this, shape, gridtools::make_gt_integer_sequence<unsigned int, 6>());
      case 7:
        return storage_resizer::apply(*this, shape, gridtools::make_gt_integer_sequence<unsigned int, 7>());
      case 8:
        return storage_resizer::apply(*this, shape, gridtools::make_gt_integer_sequence<unsigned int, 8>());
      case 9:
        return storage_resizer::apply(*this, shape, gridtools::make_gt_integer_sequence<unsigned int, 9>());
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

private:
  std::unique_ptr< DataStoreInterface>  data_store_;

public:
  ArraySpec& spec() {return spec_;}

private: // methods
  ArraySpec spec_;
};

//------------------------------------------------------------------------------

template<typename DATA_TYPE>
class ArrayT : public Array  {
public:

public:

  template<typename DataStore, typename = typename std::enable_if<gridtools::is_data_store<DataStore>::value >::type >
  ArrayT(DataStore* ds): owned_(true), data_(static_cast<void*>(ds)), Array(ds) {}

  virtual array::DataType datatype() const { return array::DataType::create<DATA_TYPE>(); }

  size_t sizeof_data() const {return sizeof(DATA_TYPE);}
  virtual void* storage() { return data_;}
  virtual const void* storage() const { return data_;}

  virtual void dump(std::ostream& os) const {
    os << "\nWARNING: Array dump not implemented\n";
  }

private:
  bool owned_;
  void* data_;
};

//------------------------------------------------------------------------------

} // namespace array
} // namespace atlas

#include "atlas/array/GridToolsArray_impl.h"
