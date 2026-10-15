#pragma once

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace atlas::array::introspection {

template <typename>
[[maybe_unused]] inline static constexpr bool always_false_v = false;


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
    if constexpr (detail::has_rank<View>::value) {
        return View::rank();
    }
    else {
        return 1;
    }
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
    using type = typename View::element_type;
};

}

template <typename View>
using element_t = typename detail::element_type<View>::type;

namespace detail {
template <typename View, typename = std::void_t<> >
struct accessor_type {
    using type = default_accessor<element_t<View>>;
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

template <typename View, typename = std::void_t<> >
struct has_static_extent : std::false_type {};

template <typename View>
struct has_static_extent<View, std::void_t<decltype(View::static_extent(std::declval<std::size_t>()))>>
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

template <typename View, typename std::enable_if_t<!has_static_extent<View>::value, int> = 0>
constexpr std::size_t get_static_extent(std::size_t /*dimension*/) noexcept {
    return dynamic_extent;  // No static extent information available, return dynamic_extent
}

template <typename View, typename std::enable_if_t<has_static_extent<View>::value, int> = 0>
constexpr std::size_t get_static_extent(std::size_t dimension) noexcept {
    return View::static_extent(dimension);  // Use mdspan's static_extent() method
}

template <typename View>
auto get_stride(const View& view, std::size_t dimension) -> decltype(view.stride(dimension)) {
    return view.stride(dimension);
}

template <std::size_t Rank>
bool valid_layout_stride(const std::array<std::size_t, Rank>& extents,
                         const std::array<std::size_t, Rank>& strides) {
    std::array<std::size_t, Rank> dimensions{};
    bool empty = false;
    for (std::size_t dimension = 0; dimension < Rank; ++dimension) {
        if (strides[dimension] == 0) {
            return false;
        }
        empty = empty || extents[dimension] == 0;
        dimensions[dimension] = dimension;
    }
    if (empty) {
        return true;
    }
    for (std::size_t position = 1; position < Rank; ++position) {
        auto dimension = dimensions[position];
        auto previous = position;
        while (previous > 0 && strides[dimensions[previous - 1]] > strides[dimension]) {
            dimensions[previous] = dimensions[previous - 1];
            --previous;
        }
        dimensions[previous] = dimension;
    }
    std::size_t previous_stride = 1;
    std::size_t previous_extent = 1;
    for (auto dimension : dimensions) {
        if (extents[dimension] > 1) {
            if (strides[dimension] / previous_stride < previous_extent) {
                return false;
            }
            previous_stride = strides[dimension];
            previous_extent = extents[dimension];
        }
    }
    return true;
}

template <typename View, std::size_t... Dimensions>
bool has_valid_layout_stride(const View& view, std::index_sequence<Dimensions...>) {
    if (!((get_stride(view, Dimensions) > 0) && ...)) {
        return false;
    }
    return valid_layout_stride(
        std::array<std::size_t, sizeof...(Dimensions)>{static_cast<std::size_t>(get_extent(view, Dimensions))...},
        std::array<std::size_t, sizeof...(Dimensions)>{static_cast<std::size_t>(get_stride(view, Dimensions))...});
}

template <typename Mapping, std::size_t Rank>
bool valid_contiguous_layout(const Mapping& mapping,
                             const std::array<std::size_t, Rank>& strides) {
    static_assert(Mapping::extents_type::rank() == Rank);
    for (std::size_t dimension = 0; dimension < Rank; ++dimension) {
        if (mapping.extents().extent(dimension) == 0) {
            return true;
        }
    }
    for (std::size_t dimension = 0; dimension < Rank; ++dimension) {
        if (mapping.extents().extent(dimension) > 1 &&
            strides[dimension] != static_cast<std::size_t>(mapping.stride(dimension))) {
            return false;
        }
    }
    return true;
}

template <typename Layout, typename View, std::size_t... I>
bool has_contiguous_layout_strides(const View& view, std::index_sequence<I...>) {
    constexpr std::size_t Rank = rank<View>();
    using Extents = dextents<std::size_t, Rank>;
    const Extents extents{static_cast<std::size_t>(get_extent(view, I))...};
    const std::array<std::size_t, Rank> strides{static_cast<std::size_t>(get_stride(view, I))...};
    const typename Layout::template mapping<Extents> mapping{extents};
    return valid_contiguous_layout(mapping, strides);
}

}  // namespace atlas::array::introspection
