#pragma once

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace atlas {
namespace mdspan_introspection_detail {

template <typename View>
using view_type_t = std::remove_reference_t<View>;

template <typename View, typename = std::void_t<> >
struct has_view_rank_member : std::false_type {};

template <typename View>
struct has_view_rank_member<View, std::void_t<decltype(view_type_t<View>::RANK)>> : std::true_type {};

template <typename View, typename = std::void_t<> >
struct has_view_rank_function : std::false_type {};

template <typename View>
struct has_view_rank_function<View, std::void_t<decltype(view_type_t<View>::rank())>> : std::true_type {};

template <typename View>
static constexpr bool has_view_rank_v = has_view_rank_member<View>::value || has_view_rank_function<View>::value;

template <typename View, typename = std::void_t<>>
struct view_rank : std::integral_constant<std::size_t, view_type_t<View>::RANK> {};

template <typename View>
struct view_rank<View, std::void_t<decltype(view_type_t<View>::rank())>>
    : std::integral_constant<std::size_t, view_type_t<View>::rank()> {};

template <typename View>
static constexpr std::size_t view_rank_v = view_rank<View>::value;

template <typename View, typename = std::void_t<> >
struct view_extents {
    using type = dextents<std::size_t, view_rank_v<View>>;
};

template <typename View>
struct view_extents<View, std::void_t<typename view_type_t<View>::extents_type>> {
    using type = typename view_type_t<View>::extents_type;
};

template <typename View>
using view_extents_t = typename view_extents<View>::type;

template <typename View, typename = std::void_t<>>
struct view_value {
    using type = std::remove_pointer_t<std::remove_reference_t<decltype(std::declval<View&>().data())>>;
};

template <typename View>
struct view_value<View, std::void_t<decltype(std::declval<View&>().data_handle())>> {
    using type = std::remove_pointer_t<std::remove_reference_t<decltype(std::declval<View&>().data_handle())>>;
};

template <typename View>
using view_value_t = typename view_value<View>::type;

template <typename View, typename = std::void_t<> >
struct view_accessor {
    using type = restrict_accessor<view_value_t<View>>;
};

template <typename View>
struct view_accessor<View, std::void_t<typename view_type_t<View>::accessor_type>> {
    using type = typename view_type_t<View>::accessor_type;
};

template <typename View>
using view_accessor_t = typename view_accessor<View>::type;

template <typename View, typename = std::void_t<> >
struct view_layout {
    using type = layout_stride;
};

template <typename View>
struct view_layout<View, std::void_t<typename view_type_t<View>::layout_type>> {
    using type = typename view_type_t<View>::layout_type;
};

template <typename View>
using view_layout_t = typename view_layout<View>::type;

template <typename>
[[maybe_unused]] inline static constexpr bool always_false_v = false;

template <typename View>
auto data_handle(View& view) -> decltype(view.data_handle()) {
    return view.data_handle();
}

template <typename View>
auto data_handle(View& view) -> decltype(view.data()) {
    return view.data();
}

template <typename View, typename = std::void_t<>>
struct has_shape : std::false_type {};

template <typename View>
struct has_shape<View, std::void_t<decltype(std::declval<const View&>().shape(std::declval<std::size_t>()))>>
    : std::true_type {};

template <typename View, typename std::enable_if_t<has_shape<View>::value, int> = 0>
auto extent(const View& view, std::size_t dimension) -> decltype(view.shape(dimension)) {
    return view.shape(dimension);
}

template <typename View, typename std::enable_if_t<!has_shape<View>::value, int> = 0>
auto extent(const View& view, std::size_t dimension) -> decltype(view.extent(dimension)) {
    return view.extent(dimension);
}

template <typename View>
auto stride(const View& view, std::size_t dimension) -> decltype(view.stride(dimension)) {
    return view.stride(dimension);
}

template <typename View, std::size_t... I>
bool has_layout_right_strides(const View& view, std::index_sequence<I...>) {
    constexpr std::size_t Rank = view_rank_v<View>;
    std::array<std::size_t, Rank> extents{static_cast<std::size_t>(extent(view, I))...};
    std::array<std::size_t, Rank> strides{static_cast<std::size_t>(stride(view, I))...};

    if (strides[Rank - 1] != 1) {
        return false;
    }

    std::size_t expected_stride{1};
    for (std::size_t i = Rank - 1; i > 0; --i) {
        expected_stride *= extents[i];
        if (strides[i - 1] != expected_stride) {
            return false;
        }
    }
    return true;
}

template <typename View, std::size_t... I>
bool dimension_start_strides_are_aligned(const View& view, std::size_t alignment, std::index_sequence<I...>) {
    using Value = std::remove_cv_t<view_value_t<View>>;
    return (((static_cast<std::size_t>(stride(view, I)) * sizeof(Value)) % alignment == 0) && ...);
}

template <std::size_t Dim, typename View, std::size_t... I>
bool dimension_start_strides_are_aligned(const View& view, std::size_t alignment, std::index_sequence<I...>) {
    return dimension_start_strides_are_aligned(
        view, alignment, std::index_sequence<(I < Dim ? I : I + 1)...>{});
}

}  // namespace mdspan_introspection_detail
}  // namespace atlas