#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "atlas/mdspan/mdspan_introspection_detail.h"

namespace atlas::array::introspection {

template <typename View>
bool can_use_layout_right(const View& view) {
    namespace introspection = array::introspection;
    constexpr std::size_t Rank = introspection::rank<View>();
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
bool can_use_layout(const View& view) {
    if constexpr (std::is_same_v<Layout, layout_right>) {
        return can_use_layout_right(view);
    }
    else if constexpr (std::is_same_v<Layout, layout_stride>) {
        return true;
    }
    else {
        static_assert(always_false_v<Layout>, "can_use_layout() is only implemented for layout_right and layout_stride");
        return false;
    }
}

template <typename View>
bool is_aligned(View& view, std::size_t alignment) {
    namespace introspection = array::introspection;
    if (alignment == 0) {
        return false;
    }
    auto address = reinterpret_cast<std::uintptr_t>(introspection::data_handle(view));
    return address % alignment == 0;
}

template <std::size_t ByteAlignment, typename View>
bool is_aligned(View& view) {
    return is_aligned(view, ByteAlignment);
}

template <std::size_t Dim, typename View>
bool is_dimension_aligned(View& view, std::size_t alignment) {
    namespace introspection = array::introspection;
    constexpr std::size_t Rank = introspection::rank<View>();

    static_assert(Dim < Rank, "Dimension must be smaller than view rank.");

    if (!is_aligned(view, alignment)) {
        return false;
    }
    if constexpr (Rank <= 1) {
        return true;
    }
    else {
        return introspection::dimension_start_strides_are_aligned<Dim>(
            view, alignment, std::make_index_sequence<Rank - 1>{});
    }
}

template <typename View>
bool is_last_dimension_aligned(View& view, std::size_t alignment) {
    constexpr std::size_t Rank = rank<View>();
    return is_dimension_aligned<Rank - 1>(view, alignment);
}


}  // namespace atlas