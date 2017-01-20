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
#include <array>
#include <type_traits>
#include "atlas/array/array_fwd.h"
#include "atlas/array/Array.h"
#include "atlas/array/ArrayUtil.h"
#include "atlas/array/DataType.h"
#include "atlas/array/gridtools/GridToolsTraits.h"
#include "eckit/exception/Exceptions.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace gridtools {

template <unsigned int Rank>
std::array<unsigned int, Rank> get_array_from_vector(std::vector<size_t> const& values) {
    std::array<unsigned int, Rank> array;
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

template<unsigned int Rank>
struct default_layout_t {

    template<typename T>
    struct get_layout;

    template<typename UInt, UInt ... Indices>
    struct get_layout<::gridtools::gt_integer_sequence<UInt, Indices...> >
    {
        using type = ::gridtools::layout_map<Indices...>;
    };

    using type = typename get_layout< typename ::gridtools::make_gt_integer_sequence<unsigned int, Rank>::type >::type;
};


  template <typename Value, typename LayoutMap>
  struct get_layout_map_component {
    // TODO: static_assert( ::gridtools::is_layout_map<LayoutMap>(), "Error: not a layout_map" );
    template <int Idx>
    struct get_component {
      GT_FUNCTION
      constexpr get_component() {}

      GT_FUNCTION constexpr static Value apply() {
        return LayoutMap::template at<Idx>();
      }
    };
  };

  template <typename Value, typename RankS>
  struct get_stride_component {
    template <int Idx>
    struct get_component {
      GT_FUNCTION
      constexpr get_component() {}

      template <typename StorageInfoPtr>
      GT_FUNCTION constexpr static Value apply(StorageInfoPtr a) {
        static_assert((::gridtools::is_storage_info<typename std::remove_pointer<StorageInfoPtr>::type>::value),
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
      GT_FUNCTION constexpr static size_t apply(StorageInfoPtr a) {
          static_assert((::gridtools::is_storage_info<typename std::remove_pointer<StorageInfoPtr>::type >::value ), "Error: not a storage_info");
          return a->template dim<Idx>();
      }
  };

  //indirection around C++11 sizeof... since it is buggy for nvcc and cray
  template<typename ...T>
  struct get_pack_size {
      using type = ::gridtools::static_uint< sizeof...(T) >;
  };

  template <typename Value, typename LayoutMap, typename... UInts>
  static gridtools::storage_traits::data_store_t<
      Value,
      gridtools::storage_traits::storage_info_t<
          0,
          get_pack_size<UInts...>::type::value,
          typename ::gridtools::zero_halo<get_pack_size<UInts...>::type::value>::type,
          LayoutMap
      >
   >*
   create_gt_storage(UInts... dims) {
      static_assert((sizeof...(dims) > 0), "Error: can not create storages without any dimension");

      constexpr static unsigned int rank = get_pack_size<UInts...>::type::value;
      typedef gridtools::storage_traits::storage_info_t<
          0,
          rank,
          typename ::gridtools::zero_halo<rank>::type,
          LayoutMap
      > storage_info_ty;
      typedef gridtools::storage_traits::data_store_t<Value, storage_info_ty> data_store_t;

      storage_info_ty si(dims...);
      data_store_t* ds = new data_store_t(si);
      ds->allocate();
      return ds;
  }

  template <typename Value, unsigned int Rank>
  static gridtools::storage_traits::data_store_t<
      Value, gridtools::storage_traits::storage_info_t<
                 0, Rank,
                 typename ::gridtools::zero_halo<Rank>::type,
                 typename default_layout_t<Rank>::type > >*
  wrap_gt_storage(
      Value* data,
      std::array<unsigned int, Rank>&& shape, std::array<unsigned int, Rank>&& strides)
  {
      static_assert((Rank > 0), "Error: can not create storages without any dimension");
      typedef gridtools::storage_traits::storage_info_t<
          0, Rank, typename ::gridtools::zero_halo<Rank>::type,
          typename default_layout_t<Rank>::type> storage_info_ty;
      typedef gridtools::storage_traits::data_store_t<Value, storage_info_ty> data_store_t;

      storage_info_ty si(shape, strides);
      data_store_t* ds = new data_store_t(si, data);
      return ds;
  }


  template<typename DataStore, typename ... Dims>
  ArraySpec make_spec(DataStore* gt_data_store_ptr, Dims...dims) {
      static_assert((::gridtools::is_data_store<DataStore>::value), "Internal Error: passing a non GT data store");

      auto storage_info_ptr = gt_data_store_ptr->get_storage_info_ptr();
      using Layout = typename DataStore::storage_info_t::Layout;

      using seq =
          ::gridtools::apply_gt_integer_sequence<typename ::gridtools::make_gt_integer_sequence<int, sizeof...(dims)>::type>;

      return ArraySpec(
          ArrayShape{(unsigned long)dims...},
          seq::template apply<
                        ArrayStrides,
                        get_stride_component<unsigned long, typename get_pack_size<Dims...>::type>::template get_component>(
                        storage_info_ptr),
          seq::template apply<
                        ArrayLayout,
                        get_layout_map_component<unsigned long, Layout>::template get_component>()
        );
  }

  //------------------------------------------------------------------------------

  template <unsigned int Rank>
  struct array_initializer;

  template<unsigned int PartDim>
  struct array_initializer_partitioned;

  //------------------------------------------------------------------------------

  template <typename Value, unsigned int Rank, unsigned int Dim>
  struct array_initializer_impl {


    static void apply(Array const& orig, Array& array_resized) {
      array_initializer_impl<Value, Rank, Dim>::apply(
            make_view<Value, Rank>(orig),
            make_view<Value, Rank>(array_resized) );
    }


    template <typename ... DimIndex>
    static void apply(ArrayView<Value, Rank> const&& orig, ArrayView<Value, Rank>&& array_resized, DimIndex... idxs) {

      for(size_t i=0; i < orig.shape(Dim); ++i) {
        array_initializer_impl<Value, Rank, Dim+1>::apply(
              std::move(orig),
              std::move(array_resized),
              idxs...,
              i );
      }

    }


  };

  //------------------------------------------------------------------------------

  template <typename Value, unsigned int Rank>
  struct array_initializer_impl<Value, Rank, Rank> {

    template <typename ... DimIndex>
    static void apply(ArrayView<Value, Rank> const&& orig,ArrayView<Value, Rank>&& array_resized, DimIndex... idxs) {
        array_resized(idxs...) = orig(idxs...);
    }

  };

  //------------------------------------------------------------------------------

  template <unsigned int Rank>
  struct array_initializer {
    template <typename... DimIndex>
    static void apply(Array const& orig, Array& array_resized, DimIndex... idxs) {
      switch (orig.datatype().kind()) {
        case DataType::KIND_REAL64: return array_initializer_impl<double,        Rank, 0>::apply(orig, array_resized, idxs...);
        case DataType::KIND_REAL32: return array_initializer_impl<float,         Rank, 0>::apply(orig, array_resized, idxs...);
        case DataType::KIND_INT32:  return array_initializer_impl<int,           Rank, 0>::apply(orig, array_resized, idxs...);
        case DataType::KIND_INT64:  return array_initializer_impl<long,          Rank, 0>::apply(orig, array_resized, idxs...);
        case DataType::KIND_UINT64: return array_initializer_impl<unsigned long, Rank, 0>::apply(orig, array_resized, idxs...);
        default: {
          std::stringstream err;
          err << "data kind " << orig.datatype().kind() << " not recognised.";
          throw eckit::BadParameter(err.str(), Here());
        }
      }
    }
  };

  //------------------------------------------------------------------------------

  template<typename Value, unsigned int Rank, unsigned int Dim, unsigned int PartDim>
  struct array_initializer_partitioned_val_impl {
    static void apply(Array const& orig, Array& dest, unsigned int pos, unsigned int offset) {
        array_initializer_partitioned_val_impl<Value, Rank, Dim, PartDim>::apply(make_view<Value, Rank>(orig), make_view<Value, Rank>(dest), pos, offset);
    }
    template <typename ... DimIndex>
    static void apply(ArrayView<Value, Rank> const&& orig, ArrayView<Value, Rank>&& dest, unsigned int pos, unsigned int offset, DimIndex... idxs) {
        for(size_t i=0; i < orig.shape(Dim); ++i)
        {
            unsigned int displ = i;
            if(Dim == PartDim && i >= pos) {
                displ += offset;
            }
            std::pair<int,int> pair_idx{i,displ};
            array_initializer_partitioned_val_impl<Value, Rank, Dim+1, PartDim>::apply(std::move(orig), std::move(dest), pos, offset, idxs..., pair_idx);
        }
    }

  };

  //------------------------------------------------------------------------------

  template <typename Value, unsigned int Rank, unsigned int PartDim>
  struct array_initializer_partitioned_val_impl<Value, Rank, Rank, PartDim> {
      template <typename ... DimIndex>
    static void apply(ArrayView<Value, Rank> const&& orig, ArrayView<Value, Rank>&& dest, unsigned int pos, unsigned int offset, DimIndex... idxs) {
        dest(std::get<1>(idxs)...) = orig(std::get<0>(idxs)...);
    }
  };

  //------------------------------------------------------------------------------

  template<unsigned int Rank, unsigned int PartDim>
  struct array_initializer_partitioned_impl {
    static void apply( Array const& orig, Array& dest, unsigned int pos, unsigned int offset) {
        switch (orig.datatype().kind()) {
          case DataType::KIND_REAL64: return array_initializer_partitioned_val_impl<double, Rank, 0, PartDim>::apply(orig, dest, pos, offset);
          case DataType::KIND_REAL32: return array_initializer_partitioned_val_impl<float, Rank, 0, PartDim>::apply(orig, dest,pos, offset );
          case DataType::KIND_INT32:  return array_initializer_partitioned_val_impl<int, Rank, 0, PartDim>::apply(orig, dest, pos, offset);
          case DataType::KIND_INT64:  return array_initializer_partitioned_val_impl<long, Rank, 0, PartDim>::apply(orig, dest, pos, offset);
          case DataType::KIND_UINT64: return array_initializer_partitioned_val_impl<unsigned long, Rank, 0, PartDim>::apply(orig, dest, pos, offset);
          default: {
            std::stringstream err;
            err << "data kind " << orig.datatype().kind() << " not recognised.";
            throw eckit::BadParameter(err.str(), Here());
          }
        }
    }
  };

  //------------------------------------------------------------------------------

  template<unsigned int PartDim>
  struct array_initializer_partitioned {
    static void apply(Array const& orig, Array& dest, unsigned int pos, unsigned int offset) {
      switch (orig.rank()) {
        case 1: return array_initializer_partitioned_impl<1, PartDim>::apply(orig, dest, pos, offset);
        case 2: return array_initializer_partitioned_impl<2, PartDim>::apply(orig, dest, pos, offset);
        case 3: return array_initializer_partitioned_impl<3, PartDim>::apply(orig, dest, pos, offset);
        case 4: return array_initializer_partitioned_impl<4, PartDim>::apply(orig, dest, pos, offset);
        case 5: return array_initializer_partitioned_impl<5, PartDim>::apply(orig, dest, pos, offset);
        case 6: return array_initializer_partitioned_impl<6, PartDim>::apply(orig, dest, pos, offset);
        case 7: return array_initializer_partitioned_impl<7, PartDim>::apply(orig, dest, pos, offset);
        case 8: return array_initializer_partitioned_impl<8, PartDim>::apply(orig, dest, pos, offset);
        case 9: return array_initializer_partitioned_impl<9, PartDim>::apply(orig, dest, pos, offset);
        default: {
          std::stringstream err;
          err << "too high Rank";
          throw eckit::BadParameter(err.str(), Here());
        }
      }
    }
  };

  //------------------------------------------------------------------------------



} // namespace gridtools
} // namespace array
} // namespace atlas
