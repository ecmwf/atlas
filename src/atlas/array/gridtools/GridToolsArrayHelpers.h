/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/array/ArrayDataStore.h"
#include "atlas/array/DataType.h"
#include "atlas/array/gridtools/GridToolsTraits.h"
#include "atlas/array_fwd.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace gridtools {

template <unsigned int Rank>
std::array<idx_t, Rank> get_array_from_vector(std::vector<idx_t> const& values) {
    std::array<idx_t, Rank> array;
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
            throw_Exception(err.str(), Here());
        }
        check_dimension_lengths_impl<TotalDims, Dim + 1>::apply(shape, d...);
    }
};

template <unsigned int TotalDims, unsigned int Dim>
struct check_dimension_lengths_impl<TotalDims, Dim,
                                    typename std::enable_if<(Dim == TotalDims - 1) && (TotalDims > 1)>::type> {
    template <typename FirstDim>
    static void apply(ArrayShape const& shape, FirstDim first_dim) {
        if (first_dim < shape[Dim]) {
            std::stringstream err;
            err << "Attempt to resize array with original size for dimension " << Dim - 1 << " of " << shape[Dim - 1]
                << " by " << first_dim << std::endl;
            throw_Exception(err.str(), Here());
        }
    }
};

template <unsigned int TotalDims, unsigned int Dim>
struct check_dimension_lengths_impl<TotalDims, Dim,
                                    typename std::enable_if<(Dim == TotalDims - 1) && (TotalDims == 1)>::type> {
    template <typename FirstDim>
    static void apply(ArrayShape const& shape, FirstDim first_dim) {
        if (first_dim < shape[Dim]) {
            std::stringstream err;
            err << "Attempt to resize array with original size for dimension " << Dim << " of " << shape[Dim] << " by "
                << first_dim << std::endl;
            throw_Exception(err.str(), Here());
        }
    }
};

template <typename... Dims>
void check_dimension_lengths(ArrayShape const& shape, Dims... d) {
    check_dimension_lengths_impl<sizeof...(d), 0>::apply(shape, d...);
}

template <unsigned int Rank>
struct default_layout_t {
    template <typename T>
    struct get_layout;

    template <typename UInt, UInt... Indices>
    struct get_layout<std::integer_sequence<UInt, Indices...>> {
        using type = ::gridtools::layout_map<Indices...>;
    };

    using type = typename get_layout<std::make_integer_sequence<::gridtools::uint_t, Rank>>::type;
};


template <typename Value, typename LayoutMap>
struct get_layout_map_component {
    // TODO: static_assert( ::gridtools::is_layout_map<LayoutMap>(), "Error: not a
    // layout_map" );
    template <int Idx>
    struct get_component {
        ATLAS_HOST_DEVICE
        constexpr get_component() {}

        ATLAS_HOST_DEVICE constexpr static Value apply() { return LayoutMap::template at<Idx>(); }
    };
};

template <typename Value>
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

template <typename Value>
struct get_shape_component {
    template <int Idx>
    struct get_component {
        ATLAS_HOST_DEVICE
        constexpr get_component() {}

        template <typename StorageInfoPtr>
        ATLAS_HOST_DEVICE constexpr static Value apply(StorageInfoPtr a) {
            static_assert((::gridtools::is_storage_info<typename std::remove_pointer<StorageInfoPtr>::type>::value),
                          "Error: not a storage_info");
            return a->template total_length<Idx>();
        }
    };
};

// indirection around C++11 sizeof... since it is buggy for nvcc and cray
template <typename... T>
struct get_pack_size {
    using type = ::gridtools::static_uint<sizeof...(T)>;
};

template <typename Value, typename LayoutMap, typename Alignment, size_t Rank>
struct gt_storage_t {
    using type =
        gridtools::storage_traits::data_store_t<Value, gridtools::storage_traits::custom_layout_storage_info_align_t<
                                                           0, LayoutMap, ::gridtools::zero_halo<Rank>, Alignment>>;
};

template <typename Value, typename LayoutMap, typename Alignment, typename... UInts>
typename gt_storage_t<Value, LayoutMap, Alignment, get_pack_size<UInts...>::type::value>::type* create_gt_storage(
    UInts... dims) {
    static_assert((sizeof...(dims) > 0), "Error: can not create storages without any dimension");
    constexpr static unsigned int rank = get_pack_size<UInts...>::type::value;
    typedef gridtools::storage_traits::custom_layout_storage_info_align_t<0, LayoutMap, ::gridtools::zero_halo<rank>,
                                                                          Alignment>
        storage_info_ty;
    typedef gridtools::storage_traits::data_store_t<Value, storage_info_ty> data_store_t;

    data_store_t* ds;
    if (::gridtools::accumulate(::gridtools::multiplies(), dims...) == 0) {
        ds = new data_store_t();
    }
    else {
        storage_info_ty si(dims...);
        ds = new data_store_t(si);
    }
    return ds;
}

template <typename Value, typename LayoutMap, typename... UInts>
typename gt_storage_t<Value, LayoutMap, ::gridtools::alignment<1>, get_pack_size<UInts...>::type::value>::type*
create_gt_storage(UInts... dims) {
    return create_gt_storage<Value, LayoutMap, ::gridtools::alignment<1>>(dims...);
}


template <typename Value, size_t Rank>
struct gt_wrap_storage_t {
    using type = gridtools::storage_traits::data_store_t<Value, gridtools::storage_traits::storage_info_t<0, Rank>>;
};

template <typename Value, unsigned int Rank>
static typename gt_wrap_storage_t<Value, Rank>::type* wrap_gt_storage(Value* data, std::array<idx_t, Rank>&& shape,
                                                                      std::array<idx_t, Rank>&& strides) {
    static_assert((Rank > 0), "Error: can not create storages without any dimension");
    typedef gridtools::storage_traits::storage_info_t<0, Rank, ::gridtools::zero_halo<Rank>> storage_info_ty;
    typedef gridtools::storage_traits::data_store_t<Value, storage_info_ty> data_store_t;
    ::gridtools::array<::gridtools::uint_t, Rank> _shape;
    ::gridtools::array<::gridtools::uint_t, Rank> _strides;
    for (unsigned int i = 0; i < Rank; ++i) {
        _shape[i]   = shape[i];
        _strides[i] = strides[i];
    }
    storage_info_ty si(_shape, _strides);
    data_store_t* ds = new data_store_t(si, data);

    return ds;
}

constexpr idx_t zero(idx_t) {
    return 0;
}

template <idx_t... Is>
ArrayStrides make_null_strides(std::integer_sequence<idx_t, Is...>) {
    return make_strides(zero(Is)...);
}

template <typename UInt>
struct my_apply_gt_integer_sequence {
    template <typename Container, template <UInt T> class Lambda, typename... ExtraTypes>
    ATLAS_HOST_DEVICE static constexpr Container apply(ExtraTypes const&... args_) {
        static_assert((std::is_same<Container, Container>::value),
                      "ERROR: my_apply_gt_integer_sequence only accepts a "
                      "std::integer_sequence type. Check the call");
        return Container(args_...);
    }
};

template <typename UInt, UInt... Indices>
struct my_apply_gt_integer_sequence<std::integer_sequence<UInt, Indices...>> {
    /**
         @brief duplicated interface for the case in which the container is an
     aggregator
       */
    template <typename Container, template <UInt T> class Lambda, typename... ExtraTypes>
    ATLAS_HOST_DEVICE static constexpr Container apply(ExtraTypes const&... args_) {
        return Container{Lambda<Indices>::apply(args_...)...};
    }
};

#ifndef __CUDACC__

template <typename DataStore, typename... Dims>
ArraySpec ATLAS_HOST make_spec(DataStore* gt_data_store_ptr, Dims... dims) {
    static_assert((::gridtools::is_data_store<DataStore>::value), "Internal Error: passing a non GT data store");

    if (gt_data_store_ptr->valid()) {
        auto storage_info_ptr = gt_data_store_ptr->get_storage_info_ptr().get();
        using Layout          = typename DataStore::storage_info_t::layout_t;
        using Alignment       = typename DataStore::storage_info_t::alignment_t;

        using seq = my_apply_gt_integer_sequence<std::make_integer_sequence<int, sizeof...(dims)>>;

        ArraySpec spec(
            ArrayShape{(idx_t)dims...},
            seq::template apply<ArrayStrides, get_stride_component<idx_t>::template get_component>(storage_info_ptr),
            seq::template apply<ArrayLayout, get_layout_map_component<idx_t, Layout>::template get_component>(),
            ArrayAlignment(Alignment::value));
        ATLAS_ASSERT(spec.allocatedSize() == storage_info_ptr->padded_total_length());
        return spec;
    }
    else {
        return ArraySpec(make_shape({dims...}),
                         make_null_strides(std::make_integer_sequence<idx_t, sizeof...(dims)>()));
    }
}
#endif

//------------------------------------------------------------------------------

}  // namespace gridtools
}  // namespace array
}  // namespace atlas
