/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <vector>
#include <iostream>
#include <type_traits>

#include "atlas/array/Array.h"
#include "atlas/array/gridtools/GridToolsTraits.h"
#include "atlas/array/gridtools/GridToolsDataStore.h"
#include "atlas/array/gridtools/GridToolsArrayHelpers.h"
#include "atlas/array/gridtools/GridToolsArrayResize.h"
#include "atlas/array/gridtools/GridToolsArrayInsert.h"

#include "atlas/array/MakeView.h"

//------------------------------------------------------------------------------

using namespace atlas::array::gridtools;

namespace atlas {
namespace array {

//------------------------------------------------------------------------------

template <typename Value,
          typename Layout,
          typename... UInts,
          typename = ::gridtools::all_integers<UInts...> >
static Array* create_array_with_layout(UInts... dims) {
  auto gt_data_store_ptr = create_gt_storage<Value, Layout>(dims...);
  using data_store_t = typename std::remove_pointer<decltype(gt_data_store_ptr)>::type;
  ArrayDataStore* data_store = new gridtools::GridToolsDataStore<data_store_t>(gt_data_store_ptr);
  return ArrayBackend<Value>::create(data_store,make_spec(gt_data_store_ptr, dims...));
}

template <typename Value,
          typename... UInts,
          typename = ::gridtools::all_integers<UInts...> >
static Array* create_array(UInts... dims) {
  auto gt_data_store_ptr = create_gt_storage<Value, typename default_layout_t<sizeof...(dims)>::type>(dims...);
  using data_store_t = typename std::remove_pointer<decltype(gt_data_store_ptr)>::type;
  ArrayDataStore* data_store = new gridtools::GridToolsDataStore<data_store_t>(gt_data_store_ptr);
  return ArrayBackend<Value>::create(data_store,make_spec(gt_data_store_ptr, dims...));
}

template <typename Value, template <class> class Storage, typename StorageInfo>
static Array* wrap_array(::gridtools::data_store< Storage<Value>, StorageInfo> * ds, const ArraySpec& spec) {
  assert(ds);
  using data_store_t = typename std::remove_pointer<decltype(ds)>::type;
  ArrayDataStore* data_store = new gridtools::GridToolsDataStore<data_store_t>(ds);
  return ArrayBackend<Value>::create(data_store,spec);
}

template <typename Value> Array* ArrayBackend<Value>::create(ArrayDataStore* datastore, const ArraySpec& spec) {
  return new ArrayT<Value>(datastore,spec);
}

template <typename Value> Array* ArrayBackend<Value>::create(const ArrayShape& shape) {
  return new ArrayT<Value>(shape);
}

template<typename Value> Array *ArrayBackend<Value>::create(const ArrayShape& shape, const ArrayLayout& layout) {
  return new ArrayT<Value>(shape,layout);
}

template <typename Value> Array* ArrayBackend<Value>::create(size_t dim0) {
  return create_array<Value>(dim0);
}
template <typename Value> Array* ArrayBackend<Value>::create(size_t dim0, size_t dim1) {
  return create_array<Value>(dim0,dim1);
}
template <typename Value> Array* ArrayBackend<Value>::create(size_t dim0, size_t dim1, size_t dim2) {
  return create_array<Value>(dim0,dim1,dim2);
}
template <typename Value> Array* ArrayBackend<Value>::create(size_t dim0, size_t dim1, size_t dim2, size_t dim3) {
  return create_array<Value>(dim0,dim1,dim2,dim3);
}
template <typename Value> Array* ArrayBackend<Value>::create(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4) {
  return create_array<Value>(dim0,dim1,dim2,dim3,dim4);
}

template <typename Value> Array* ArrayBackend<Value>::wrap(Value* data, const ArrayShape& shape) {
  return wrap(data,ArraySpec(shape));
}

template <typename Value> Array* ArrayBackend<Value>::wrap(Value* data, const ArraySpec& spec) {
  ArrayShape const& shape = spec.shape();
  ArrayStrides const& strides = spec.strides();

  assert(shape.size() > 0);

  switch (shape.size()) {
    case 1: return wrap_array( wrap_gt_storage<Value, 1>(data, get_array_from_vector<1>(shape), get_array_from_vector<1>(strides)), spec);
    case 2: return wrap_array( wrap_gt_storage<Value, 2>(data, get_array_from_vector<2>(shape), get_array_from_vector<2>(strides)), spec);
    case 3: return wrap_array( wrap_gt_storage<Value, 3>(data, get_array_from_vector<3>(shape), get_array_from_vector<3>(strides)), spec);
    case 4: return wrap_array( wrap_gt_storage<Value, 4>(data, get_array_from_vector<4>(shape), get_array_from_vector<4>(strides)), spec);
    case 5: return wrap_array( wrap_gt_storage<Value, 5>(data, get_array_from_vector<5>(shape), get_array_from_vector<5>(strides)), spec);
    case 6: return wrap_array( wrap_gt_storage<Value, 6>(data, get_array_from_vector<6>(shape), get_array_from_vector<6>(strides)), spec);
    case 7: return wrap_array( wrap_gt_storage<Value, 7>(data, get_array_from_vector<7>(shape), get_array_from_vector<7>(strides)), spec);
    case 8: return wrap_array( wrap_gt_storage<Value, 8>(data, get_array_from_vector<8>(shape), get_array_from_vector<8>(strides)), spec);
    case 9: return wrap_array( wrap_gt_storage<Value, 9>(data, get_array_from_vector<9>(shape), get_array_from_vector<9>(strides)), spec);
    default: {
      std::stringstream err;
      err << "shape not recognized";
      throw eckit::BadParameter(err.str(), Here());
    }
  }
}

Array* Array::create(DataType datatype, const ArrayShape& shape)
{
  switch( datatype.kind() )
  {
    case DataType::KIND_REAL64: return ArrayBackend<double>::create(shape);
    case DataType::KIND_REAL32: return ArrayBackend<float>::create(shape);
    case DataType::KIND_INT32:  return ArrayBackend<int>::create(shape);
    case DataType::KIND_INT64:  return ArrayBackend<long>::create(shape);
    case DataType::KIND_UINT64: return ArrayBackend<unsigned long>::create(shape);
    default:
    {
      std::stringstream err; err << "data kind " << datatype.kind() << " not recognised.";
      throw eckit::BadParameter(err.str(),Here());
    }
  }
  return 0;
}

Array* Array::create(DataType datatype, const ArrayShape& shape, const ArrayLayout& layout)
{
  switch( datatype.kind() )
  {
  case DataType::KIND_REAL64: return ArrayBackend<double>::create(shape,layout);
  case DataType::KIND_REAL32: return ArrayBackend<float>::create(shape,layout);
  case DataType::KIND_INT32:  return ArrayBackend<int>::create(shape,layout);
  case DataType::KIND_INT64:  return ArrayBackend<long>::create(shape,layout);
  case DataType::KIND_UINT64: return ArrayBackend<unsigned long>::create(shape,layout);

    default:
    {
      std::stringstream err; err << "data kind " << datatype.kind() << " not recognised.";
      throw eckit::BadParameter(err.str(),Here());
    }
  }
  return 0;
}

//------------------------------------------------------------------------------

template <typename Value> void ArrayT<Value>::insert(size_t idx1, size_t size1) {
  gridtools::ArrayBackendInsert(*this).insert(idx1,size1);
}

//------------------------------------------------------------------------------

template <typename Value> void ArrayT<Value>::resize(const ArrayShape& shape) {
  gridtools::ArrayBackendResize(*this).resize(shape);
}

template <typename Value> void ArrayT<Value>::resize(size_t dim0) {
  gridtools::ArrayBackendResize(*this).resize(dim0);
}

template <typename Value> void ArrayT<Value>::resize(size_t dim0, size_t dim1) {
  gridtools::ArrayBackendResize(*this).resize(dim0,dim1);
}

template <typename Value> void ArrayT<Value>::resize(size_t dim0, size_t dim1, size_t dim2) {
  gridtools::ArrayBackendResize(*this).resize(dim0,dim1,dim2);
}

template <typename Value> void ArrayT<Value>::resize(size_t dim0, size_t dim1, size_t dim2, size_t dim3) {
  gridtools::ArrayBackendResize(*this).resize(dim0,dim1,dim2,dim3);
}

template <typename Value> void ArrayT<Value>::resize(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4) {
  gridtools::ArrayBackendResize(*this).resize(dim0,dim1,dim2,dim3,dim4);
}

template <typename Value> void ArrayT<Value>::dump(std::ostream& os) const {
  os << "\nWARNING: Array dump not implemented\n";
}

//------------------------------------------------------------------------------

template <typename Value>
class ArrayTConstructor {
public:
  ArrayTConstructor(ArrayT<Value>& array) : array_(array) {}

  template <typename... UInts,
            typename = ::gridtools::all_integers<UInts...> >
  void operator()(UInts... dims) {
    construct(dims...);
  }

  void operator()(const ArrayShape& shape) {
    construct(shape);
  }

  void operator()(const ArrayShape& shape, const ArrayLayout& layout) {
      construct_with_layout(shape,layout);
  }

private:

  template <typename... UInts,
            typename = ::gridtools::all_integers<UInts...> >
  void construct(UInts... dims)
  {
    auto gt_storage = create_gt_storage<Value,typename default_layout_t<sizeof...(dims)>::type>(dims...);
    using data_store_t = typename std::remove_pointer<decltype(gt_storage)>::type;
    array_.data_store_ = std::unique_ptr<ArrayDataStore>(new gridtools::GridToolsDataStore<data_store_t>(gt_storage));
    array_.spec_ = make_spec(gt_storage,dims...);
  }

  template <typename Layout,
            typename... UInts,
            typename = ::gridtools::all_integers<UInts...> >
  void construct_with_layout(UInts... dims) {
    auto gt_data_store_ptr = create_gt_storage<Value, Layout>(dims...);
    using data_store_t = typename std::remove_pointer<decltype(gt_data_store_ptr)>::type;
    array_.data_store_ = std::unique_ptr<ArrayDataStore>( new gridtools::GridToolsDataStore<data_store_t>(gt_data_store_ptr) );
    array_.spec_ = make_spec(gt_data_store_ptr, dims...);
  }

  void construct(const ArrayShape& shape) {
    assert(shape.size() > 0);
    switch (shape.size()) {
      case 1: return construct(shape[0]);
      case 2: return construct(shape[0],shape[1]);
      case 3: return construct(shape[0],shape[1],shape[2]);
      case 4: return construct(shape[0],shape[1],shape[2],shape[3]);
      case 5: return construct(shape[0],shape[1],shape[2],shape[3],shape[4]);
      default: {
        std::stringstream err;
        err << "shape not recognized";
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }

  void construct_with_layout(const ArrayShape& shape, const ArrayLayout& layout) {
    ASSERT( shape.size() > 0 );
    ASSERT( shape.size() == layout.size() );
    switch (shape.size()) {
      case 1: return construct(shape[0]);
      case 2: {
        if( layout[0]==0 && layout[1]==1 ) return construct_with_layout<::gridtools::layout_map<0,1>>(shape[0],shape[1]);
        if( layout[0]==1 && layout[1]==0 ) return construct_with_layout<::gridtools::layout_map<1,0>>(shape[0],shape[1]);
      }
      case 3: {
        if( layout[0]==0 && layout[1]==1 && layout[2] == 2) return construct_with_layout<::gridtools::layout_map<0,1,2>>(shape[0],shape[1],shape[2]);
        if( layout[0]==2 && layout[1]==1 && layout[2] == 0) return construct_with_layout<::gridtools::layout_map<2,1,0>>(shape[0],shape[1],shape[2]);
      }
      case 4: {
        if( layout[0]==0 && layout[1]==1 && layout[2] == 2 && layout[3] == 3) return construct_with_layout<::gridtools::layout_map<0,1,2,3>>(shape[0],shape[1],shape[2],shape[3]);
        if( layout[0]==3 && layout[1]==2 && layout[2] == 1 && layout[3] == 0) return construct_with_layout<::gridtools::layout_map<3,2,1,0>>(shape[0],shape[1],shape[2],shape[3]);
      }
      case 5: {
        if( layout[0]==0 && layout[1]==1 && layout[2] == 2 && layout[3] == 3 && layout[4] == 4) return construct_with_layout<::gridtools::layout_map<0,1,2,3,4>>(shape[0],shape[1],shape[2],shape[3],shape[4]);
        if( layout[0]==4 && layout[1]==3 && layout[2] == 2 && layout[3] == 1 && layout[4] == 0) return construct_with_layout<::gridtools::layout_map<4,3,2,1,0>>(shape[0],shape[1],shape[2],shape[3],shape[4]);
      }
      default: {
        std::stringstream err;
        if( shape.size() > 5 )
          err << "shape not recognized";
        else {
          err << "Layout < ";
          for( size_t j=0; j<layout.size(); ++j) err << layout[j] << " ";
          err << "> not implemented in Atlas.";
        }
        throw eckit::BadParameter(err.str(), Here());
      }
    }
  }

private:
  ArrayT<Value>& array_;
};


template <typename Value> ArrayT<Value>::ArrayT(ArrayDataStore* ds, const ArraySpec& spec) {
  data_store_ = std::unique_ptr<ArrayDataStore>(ds);
  spec_ = spec;
}

template <typename Value> ArrayT<Value>::ArrayT(size_t dim0) {
  ArrayTConstructor<Value>(*this)(dim0);
}
template <typename Value> ArrayT<Value>::ArrayT(size_t dim0, size_t dim1) {
  ArrayTConstructor<Value>(*this)(dim0,dim1);
}
template <typename Value> ArrayT<Value>::ArrayT(size_t dim0, size_t dim1, size_t dim2) {
  ArrayTConstructor<Value>(*this)(dim0,dim1,dim2);
}
template <typename Value> ArrayT<Value>::ArrayT(size_t dim0, size_t dim1, size_t dim2, size_t dim3) {
  ArrayTConstructor<Value>(*this)(dim0,dim1,dim2,dim3);
}
template <typename Value> ArrayT<Value>::ArrayT(size_t dim0, size_t dim1, size_t dim2, size_t dim3, size_t dim4) {
  ArrayTConstructor<Value>(*this)(dim0,dim1,dim2,dim3,dim4);
}

template <typename Value> ArrayT<Value>::ArrayT(const ArrayShape& shape) {
  ArrayTConstructor<Value>(*this)(shape);
}

template <typename Value> ArrayT<Value>::ArrayT(const ArrayShape& shape, const ArrayLayout& layout) {
  ArrayTConstructor<Value>(*this)(shape,layout);
}

template <typename Value> ArrayT<Value>::ArrayT(const ArraySpec& spec) {
  if( not spec.contiguous() )     NOTIMP;
  ArrayTConstructor<Value>(*this)(spec.shape(),spec.layout());
}

//------------------------------------------------------------------------------

template class ArrayBackend<int>;
template class ArrayBackend<long>;
template class ArrayBackend<float>;
template class ArrayBackend<double>;
template class ArrayBackend<unsigned long>;

template class ArrayT<int>;
template class ArrayT<long>;
template class ArrayT<float>;
template class ArrayT<double>;
template class ArrayT<unsigned long>;

} // namespace array
} // namespace atlas
