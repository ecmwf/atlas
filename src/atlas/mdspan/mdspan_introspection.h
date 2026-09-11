/*
 * (C) Copyright 2026- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

 #pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "atlas/mdspan/mdspan_introspection_detail.h"

namespace atlas::array {
    template <typename View>
    auto data_handle(View& view) {
        return introspection::data_handle(view);
    }

    template <typename View>
    auto extent(const View& view, std::size_t dimension) {
        return introspection::get_extent(view, dimension);
    }

    template <std::size_t Dim, typename View>
    auto extent(const View& view) {
        return extent(view, Dim);
    }

    template <typename View>
    auto stride(const View& view, std::size_t dimension) {
        return introspection::get_stride(view, dimension);
    }

    template <std::size_t Dim, typename View>
    auto stride(const View& view) {
        return stride(view, Dim);
    }

    template <typename View>
    auto last_extent(const View& view) {
        return extent<introspection::rank<View>() - 1>(view);
    }

    template <typename View>
    constexpr bool is_aligned(View& view, std::size_t alignment) {
        if (alignment == 0) {
            return false;
        }
        auto address = reinterpret_cast<std::uintptr_t>(data_handle(view));
        return address % alignment == 0;
    }

    template <std::size_t ByteAlignment, typename View>
    constexpr bool is_aligned(View& view) {
        return is_aligned(view, ByteAlignment);
    }

    template <typename View>
    constexpr auto rank() {
        return introspection::rank<View>();
    }

    template <typename View>
    constexpr auto rank(const View&) {
        return rank<View>();
    }

    template <typename View>
    bool conforms_layout_right(const View& view) {
        constexpr std::size_t Rank = rank<View>();
        if constexpr (std::is_same_v<introspection::layout_t<View>, layout_right>) {
            return true;
        }
        else if constexpr (Rank <= 1) {
            return true;
        }
        else {
            return introspection::has_layout_right_strides(view, std::make_index_sequence<Rank>{});
        }
    }

    template <typename Layout, typename View>
    bool conforms_layout(const View& view) {
        if constexpr (std::is_same_v<Layout, layout_right>) {
            return conforms_layout_right(view);
        }
        else if constexpr (std::is_same_v<Layout, layout_stride>) {
            return true;
        }
        else {
            static_assert(always_false_v<Layout>, "conforms_layout() is only implemented for layout_right and layout_stride");
            return false;
        }
    }

    template <std::size_t Dim, typename View>
    constexpr bool is_dimension_aligned(View& view, std::size_t alignment) {
        namespace introspection = array::introspection;
        constexpr std::size_t Rank = rank<View>();

        static_assert(Dim < Rank, "Dimension must be smaller than view rank.");

        if (!is_aligned(view, alignment)) {
            return false;
        }
        if constexpr (Rank <= 1) {
            return true;
        }
        else {
            return array::stride<Dim-1>(view) * sizeof(introspection::element_t<View>) % alignment == 0;
        }
    }

    template <typename View>
    constexpr bool is_last_dimension_aligned(View& view, std::size_t alignment) {
        constexpr std::size_t Rank = rank<View>();
        return is_dimension_aligned<Rank - 1>(view, alignment);
    }

    template <typename View>
    static constexpr bool is_mdspan() { return introspection::is_mdspan<View>(); }

    template <typename View>
    static constexpr bool is_mdspan(const View&) { return is_mdspan<View>(); }
}
