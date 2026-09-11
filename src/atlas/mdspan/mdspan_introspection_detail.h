#pragma once

#include <array>
#include <cstddef>
#include <iostream>
#include <type_traits>
#include <utility>

#ifndef DEBUG_IS_DIMENSION_ALIGNED
#define DEBUG_IS_DIMENSION_ALIGNED 1
#endif

namespace atlas {
template <typename>
[[maybe_unused]] inline static constexpr bool always_false_v = false;
}

namespace atlas::array::introspection {

namespace detail {
template <typename T, typename = void>
struct is_mdspan : std::false_type {};

template <typename T>
struct is_mdspan<T, std::void_t<
    typename T::element_type,
    typename T::data_handle_type,
    typename T::mapping_type,
    decltype(std::declval<T>().data_handle()),
    decltype(std::declval<T>().accessor())
>> : std::true_type {};
}

template <typename View>
static constexpr bool is_mdspan() { return detail::is_mdspan<View>::value; }

namespace detail {
template <typename View, typename = std::void_t<> >
struct has_rank : std::false_type {};

template <typename View>
struct has_rank<View, std::void_t<decltype(View::rank())>> : std::true_type {};
}

template <typename View>
static constexpr bool has_rank() { return detail::has_rank<View>::value; }

template <typename View>
static constexpr std::size_t rank() {
    return detail::has_rank<View>::value ? View::rank() : 1;
}

namespace detail {
template <typename View, typename = std::void_t<> >
struct has_data : std::false_type {};
template <typename View>
struct has_data<View, std::void_t<decltype(std::declval<View&>().data())>> : std::true_type {};
}


template <typename View>
static constexpr bool has_data() { return detail::has_data<View>::value; }

namespace detail {
template <typename View, typename = std::void_t<> >
struct has_data_handle : std::false_type {};
template <typename View>
struct has_data_handle<View, std::void_t<decltype(std::declval<View&>().data_handle())>> : std::true_type {};

}

template <typename View>
static constexpr bool has_data_handle() { return detail::has_data_handle<View>::value; }


namespace detail {
template <typename View, typename = std::void_t<> >
struct extents {
    using type = dextents<std::size_t, rank<View>()>;
};

template <typename View>
struct extents<View, std::void_t<typename View::extents_type>> {
    using type = typename View::extents_type;
};
}

template <typename View>
using extents_t = typename detail::extents<View>::type;


namespace detail {
template <typename View, typename = std::void_t<>>
struct element_type_from_data_or_data_handle {
    using type = std::remove_pointer_t<std::remove_reference_t<decltype(std::declval<View&>().data())>>;
};

template <typename View>
struct element_type_from_data_or_data_handle<View, std::void_t<decltype(std::declval<View&>().data_handle())>> {
    using type = std::remove_pointer_t<std::remove_reference_t<decltype(std::declval<View&>().data_handle())>>;
};

template <typename View, typename = std::void_t<>>
struct element_type {
    using type = typename element_type_from_data_or_data_handle<View>::type;
};

template <typename View>
struct element_type<View, std::void_t<typename View::element_type>> {
    using raw_type = typename View::element_type;
    using type = std::conditional_t<std::is_const_v<std::remove_reference_t<View>>,
                                     std::add_const_t<raw_type>, raw_type>;
};

}

template <typename View>
using element_t = typename detail::element_type<View>::type;

namespace detail {
template <typename View, typename = std::void_t<> >
struct accessor_type {
    using type = restrict_accessor<element_t<View>>;
};

template <typename View>
struct accessor_type<View, std::void_t<typename View::accessor_type>> {
    using type = typename View::accessor_type;
};

}
template <typename View>
using accessor_t = typename detail::accessor_type<View>::type;

namespace detail {
template <typename View, typename = std::void_t<> >
struct layout_type {
    using type = layout_stride;
};

template <typename View>
struct layout_type<View, std::void_t<typename View::layout_type>> {
    using type = typename View::layout_type;
};

}
template <typename View>
using layout_t = typename detail::layout_type<View>::type;

template <typename View, typename = std::void_t<> >
struct has_extent : std::false_type {};

template <typename View>
struct has_extent<View, std::void_t<decltype(std::declval<const View&>().extent(std::declval<std::size_t>()))>>
    : std::true_type {};

template <typename View, typename std::enable_if_t<has_data_handle<View>(), int> = 0>
auto data_handle(View& view) -> decltype(view.data_handle()) {
    return view.data_handle();
}

template <typename View, typename std::enable_if_t<!has_data_handle<View>(), int> = 0>
auto data_handle(View& view) -> decltype(view.data()) {
    return view.data();
}

template <typename View, typename std::enable_if_t<has_extent<View>::value, int> = 0>
auto get_extent(const View& view, std::size_t dimension) -> decltype(view.extent(dimension)) {
    return view.extent(dimension);
}

template <typename View>
auto get_stride(const View& view, std::size_t dimension) -> decltype(view.stride(dimension)) {
    return view.stride(dimension);
}

template <typename View, std::size_t... I>
bool has_layout_right_strides(const View& view, std::index_sequence<I...>) {
    constexpr std::size_t Rank = rank<View>();
    std::array<std::size_t, Rank> extents{static_cast<std::size_t>(get_extent(view, I))...};
    std::array<std::size_t, Rank> strides{static_cast<std::size_t>(get_stride(view, I))...};

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

}  // namespace atlas::array::introspection
