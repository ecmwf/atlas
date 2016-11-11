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
#include "atlas/array/GridToolsTraits.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

//template<unsigned int NDims>
//struct default_layout {

//    template<typename T>
//    struct get_layout;

//    template<typename UInt, UInt ... Indices>
//    struct get_layout<gridtools::gt_integer_sequence<UInt, Indices...> >
//    {
//        using type = gridtools::layout_map<Indices...>;
//    };

//    using type = typename get_layout< typename gridtools::make_gt_integer_sequence<unsigned int, NDims>::type >::type;
//};

//template <unsigned int NDims>
//std::array<unsigned int, NDims> get_array_from_vector(std::vector<size_t> const& values) {
//    std::array<unsigned int, NDims> array;
//    std::copy(values.begin(), values.end(), array.begin());
//    return array;
//}

//template <unsigned int TotalDims, unsigned int Dim, typename = void>
//struct check_dimension_lengths_impl {
//  template <typename FirstDim, typename... Dims>
//  static void apply(ArrayShape const& shape, FirstDim first_dim, Dims... d) {
//    if (first_dim < shape[Dim]) {
//      std::stringstream err;
//      err << "Attempt to resize array with original size for dimension " << Dim << " of " << shape[Dim] << " by "
//          << first_dim << std::endl;
//      throw eckit::BadParameter(err.str(), Here());
//    }
//    check_dimension_lengths_impl<TotalDims, Dim + 1>::apply(shape, d...);
//  }
//};

//template <unsigned int TotalDims, unsigned int Dim>
//struct check_dimension_lengths_impl<TotalDims, Dim, typename std::enable_if< (Dim == TotalDims-1)>::type > {
//  template <typename FirstDim>
//  static void apply(ArrayShape const& shape, FirstDim first_dim) {
//    if (first_dim < shape[Dim]) {
//      std::stringstream err;
//      err << "Attempt to resize array with original size for dimension " << Dim - 1 << " of "
//          << shape[Dim - 1] << " by " << first_dim << std::endl;
//      throw eckit::BadParameter(err.str(), Here());
//    }
//  }
//};

//template<typename ... Dims>
//void check_dimension_lengths(ArrayShape const&  shape, Dims...d) {
//    check_dimension_lengths_impl<sizeof...(d), 0>::apply(shape, d...);
//}

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

} // namespace array
} // namespace atlas
