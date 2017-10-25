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
#include "atlas/library/config.h"
#include "atlas/array_fwd.h"
#include "atlas/array.h"
#include "atlas/array/ArrayUtil.h"
#include "atlas/array/DataType.h"
#include "atlas/array/gridtools/GridToolsTraits.h"
#include "atlas/array/helpers/ArrayInitializer.h"
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
struct check_dimension_lengths_impl<TotalDims, Dim, typename std::enable_if< (Dim == TotalDims-1) && (TotalDims > 1) >::type > {
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

template <unsigned int TotalDims, unsigned int Dim>
struct check_dimension_lengths_impl<TotalDims, Dim, typename std::enable_if<
    (Dim == TotalDims-1) &&
    (TotalDims == 1) >::type > {
  template <typename FirstDim>
  static void apply(ArrayShape const& shape, FirstDim first_dim) {
    if (first_dim < shape[Dim]) {
      std::stringstream err;
      err << "Attempt to resize array with original size for dimension " << Dim << " of " << shape[Dim] << " by "
          << first_dim << std::endl;
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
      ATLAS_HOST_DEVICE
      constexpr get_component() {}

      ATLAS_HOST_DEVICE constexpr static Value apply() {
        return LayoutMap::template at<Idx>();
      }
    };
  };

  template <typename Value, typename RankS>
  struct get_stride_component {
    template <int Idx>
    struct get_component {
      ATLAS_HOST_DEVICE
      constexpr get_component() {}

      template <typename StorageInfoPtr>
      ATLAS_HOST_DEVICE constexpr static Value apply(StorageInfoPtr a) {
        static_assert((::gridtools::is_storage_info<typename std::remove_pointer<StorageInfoPtr>::type>::value),
                      "Error: not a storage_info");
        return a->template stride<Idx>();
      }
    };
  };

  template < int Idx >
  struct get_shape_component {

      ATLAS_HOST_DEVICE
      constexpr get_shape_component() {}

      template < typename StorageInfoPtr>
      ATLAS_HOST_DEVICE constexpr static size_t apply(StorageInfoPtr a) {
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


    template < typename UInt >
    struct my_apply_gt_integer_sequence {
        template < typename Container, template < UInt T > class Lambda, typename... ExtraTypes >
        ATLAS_HOST_DEVICE static constexpr Container apply(ExtraTypes const &... args_) {
            static_assert((boost::is_same< Container, Container >::value),
                "ERROR: apply_gt_integer_sequence only accepts a gt_integer_sequence type. Check the call");
            return Container(args_...);
        }
    };

    template < typename UInt, UInt... Indices >
    struct my_apply_gt_integer_sequence< ::gridtools::gt_integer_sequence< UInt, Indices... > > {

        /**
           @brief duplicated interface for the case in which the container is an aggregator
         */
        template < typename Container,
            template < UInt T > class Lambda,
            typename... ExtraTypes >
        ATLAS_HOST_DEVICE static constexpr Container apply(ExtraTypes const &... args_) {
            return Container{Lambda< Indices >::apply(args_...)...};
        }
     };



  template<typename DataStore, typename ... Dims>
  ArraySpec make_spec(DataStore* gt_data_store_ptr, Dims...dims) {
      static_assert((::gridtools::is_data_store<DataStore>::value), "Internal Error: passing a non GT data store");

      auto storage_info_ptr = gt_data_store_ptr->get_storage_info_ptr();
      using Layout = typename DataStore::storage_info_t::Layout;

      using seq = my_apply_gt_integer_sequence<typename ::gridtools::make_gt_integer_sequence<int, sizeof...(dims)>::type>;

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

} // namespace gridtools
} // namespace array
} // namespace atlas
