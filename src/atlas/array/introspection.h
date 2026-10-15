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
#include <numeric>
#include <type_traits>
#include <utility>

#include "atlas/array/mdspan/mdspan.h"
#include "atlas/array/detail/introspection_detail.h"

namespace atlas::array {

    template <typename View>
    constexpr auto rank() {
        return introspection::rank<View>();
    }

    template <typename View>
    constexpr auto rank(const View&) {
        return rank<View>();
    }

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

    // Returns either the static extent or 'dynamic_extent' (which is not the actual extent).
    template <typename View>
    constexpr auto static_extent(std::size_t dimension) {
        return introspection::get_static_extent<View>(dimension);
    }

    template <typename View>
    constexpr auto static_extent(const View&, std::size_t dimension) {
        return static_extent<View>(dimension);
    }

    template <std::size_t Dim, typename View>
    constexpr auto static_extent(const View&) {
        return static_extent<View>(Dim);
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
        return extent<rank<View>() - 1>(view);
    }

    /// Check only the view's starting address, not its strides or contiguity.
    /// Use for one specific view or subview passed to an aligned kernel.
    /// An alignment of zero is invalid and returns false.
    template <typename View>
    constexpr bool is_aligned(const View& view, std::size_t alignment) {
        if (alignment == 0) {
            return false;
        }
        auto address = reinterpret_cast<std::uintptr_t>(data_handle(view));
        return address % alignment == 0;
    }

    template <std::size_t ByteAlignment, typename View>
    constexpr bool is_aligned(const View& view) {
        return is_aligned(view, ByteAlignment);
    }

    template <typename View>
    bool conforms_layout_right(const View& view) {
        constexpr std::size_t Rank = rank<View>();
        if constexpr (std::is_same_v<introspection::layout_t<View>, layout_right>) {
            return true;
        }
        else if constexpr (std::is_same_v<introspection::layout_t<View>, layout_left>) {
            return false;
        }
        else if constexpr (Rank == 0 || !introspection::has_rank<View>()) {
            return true;
        }
        else {
            return introspection::has_contiguous_layout_strides<layout_right>(view, std::make_index_sequence<Rank>{});
        }
    }

    template <typename View>
    bool conforms_layout_left(const View& view) {
        constexpr std::size_t Rank = rank<View>();
        if constexpr (std::is_same_v<introspection::layout_t<View>, layout_left>) {
            return true;
        }
        else if constexpr (std::is_same_v<introspection::layout_t<View>, layout_right>) {
            return false;
        }
        else if constexpr (Rank == 0 || !introspection::has_rank<View>()) {
            return true;
        }
        else {
            return introspection::has_contiguous_layout_strides<layout_left>(view, std::make_index_sequence<Rank>{});
        }
    }

    template <typename Layout, typename View>
    bool conforms_layout(const View& view) {
        if constexpr (std::is_same_v<Layout, layout_right>) {
            return conforms_layout_right(view);
        }
        else if constexpr (std::is_same_v<Layout, layout_left>) {
            return conforms_layout_left(view);
        }
        else if constexpr (std::is_same_v<Layout, layout_stride>) {
            if constexpr (!introspection::has_rank<View>() || rank<View>() == 0) {
                return true;
            }
            else {
                return introspection::has_valid_layout_stride(view, std::make_index_sequence<rank<View>()>{});
            }
        }
        else {
            static_assert(introspection::always_false_v<Layout>, "conforms_layout() is only implemented for layout_left, layout_right and layout_stride");
            return false;
        }
    }

    /// Check whether advancing one index along Dim preserves byte alignment.
    /// This checks spacing only, even for singleton dimensions; the base may be unaligned.
    /// For a 2D double view with strides {8, 1}, stride 0 preserves 64-byte alignment.
    /// Combine with is_aligned() to establish aligned row starts in that example.
    /// An alignment of zero is invalid and returns false.
    template <std::size_t Dim, typename View>
    constexpr bool is_stride_aligned(const View& view, std::size_t alignment) {
        static_assert(Dim < rank<View>(), "Dimension must be smaller than view rank.");
        if (alignment == 0) {
            return false;
        }
        const auto element_period = alignment / std::gcd(alignment, sizeof(introspection::element_t<View>));
        return array::stride<Dim>(view) % element_period == 0;
    }

    namespace introspection {
        template <typename View, std::size_t... Dimensions>
        constexpr bool preceding_strides_aligned([[maybe_unused]] const View& view,
                                                 [[maybe_unused]] std::size_t alignment,
                                                 std::index_sequence<Dimensions...>) {
            return ((array::extent<Dimensions>(view) <= 1 ||
                     array::is_stride_aligned<Dimensions>(view, alignment)) && ...);
        }
    }

    /// Check alignment of every trailing subview starting at dimension Dim.
    /// Fixing any valid indices before Dim must yield an aligned starting address.
    /// Checks the base and every preceding stride whose extent exceeds one.
    /// Dim == 0 checks only the base, equivalent to is_aligned(view, alignment).
    /// Use for kernels applied to all trailing slices, including arbitrary strided views.
    /// Does not guarantee contiguity or alignment of each element inside a slice.
    /// For empty views, still checks the base and preceding non-singleton strides.
    template <std::size_t Dim, typename View>
    constexpr bool is_dimension_aligned(const View& view, std::size_t alignment) {
        static_assert(Dim < rank<View>(), "Dimension must be smaller than view rank.");
        return is_aligned(view, alignment) &&
               introspection::preceding_strides_aligned(view, alignment, std::make_index_sequence<Dim>{});
    }

    /// Check all last-dimension slice starts; pair with a layout/contiguity check
    /// before using contiguous aligned SIMD loads on each slice.
    template <typename View>
    constexpr bool is_last_dimension_aligned(const View& view, std::size_t alignment) {
        constexpr std::size_t Rank = rank<View>();
        static_assert(Rank > 0, "Last-dimension alignment requires a nonzero view rank.");
        return is_dimension_aligned<Rank - 1>(view, alignment);
    }

    template <typename View>
    static constexpr bool is_mdspan() { return introspection::is_mdspan<View>(); }

    template <typename View>
    static constexpr bool is_mdspan(const View&) { return is_mdspan<View>(); }
}
