#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "atlas/mdspan/mdspan.h"
#include "atlas/mdspan/mdspan_introspection.h"

namespace atlas {
namespace make_mdspan_helpers {

template <typename View>
static constexpr std::size_t view_rank_v = mdspan_introspection_detail::view_rank_v<View>;

template <typename View>
using view_value_t = typename mdspan_introspection_detail::view_value_t<View>;

template <typename View>
using view_extents_t = typename mdspan_introspection_detail::view_extents_t<View>;

template <typename View>
using view_layout_t = typename mdspan_introspection_detail::view_layout_t<View>;

template <typename View>
using view_accessor_t = typename mdspan_introspection_detail::view_accessor_t<View>;


template <typename Container, typename = std::void_t<> >
struct is_container_like : std::false_type {};

template <typename Container>
struct is_container_like<Container, std::void_t<
    decltype(std::declval<Container&>().data()),
    decltype(std::declval<Container&>().size())>>
    : std::integral_constant<bool, !mdspan_introspection_detail::has_view_rank_v<Container>> {};

template <typename Container>
static constexpr bool is_container_like_v = is_container_like<Container>::value;

template <typename T, std::size_t Rank>
struct extract_extents { using type = T; };

template <typename T, std::size_t N, std::size_t Rank>
struct extract_extents<std::array<T, N>, Rank> {
    using type = dextents<T, Rank>;
};

template <typename IndexType, std::size_t... ExtentsPack, std::size_t Rank>
struct extract_extents<extents<IndexType, ExtentsPack...>, Rank> {
    using type = extents<IndexType, ExtentsPack...>;
};

template <typename T, std::size_t Rank>
using extract_extents_t = typename extract_extents<T, Rank>::type;

template <typename T, typename = std::void_t<>>
struct is_extent_like : std::false_type {};

template <typename T>
struct is_extent_like<T, std::void_t<decltype(std::declval<T>().rank())>> : std::true_type {};

template <typename T>
static constexpr bool is_extent_like_v = is_extent_like<T>::value;

template <typename T>
struct is_static_extent_spec : std::false_type {};

template <std::size_t Dim, std::size_t N>
struct is_static_extent_spec<extent_with_static_dim<Dim, N>> : std::true_type {};

template <std::size_t N>
struct is_static_extent_spec<extent_with_static_last_dim<N>> : std::true_type {};

template <typename T>
static constexpr bool is_static_extent_spec_v = is_static_extent_spec<T>::value;

template <typename T, typename = std::void_t<> >
struct extents_rank;

template <typename T, std::size_t N>
struct extents_rank<std::array<T, N>> : std::integral_constant<std::size_t, N> {};

template <typename IndexType, std::size_t... ExtentsPack>
struct extents_rank<extents<IndexType, ExtentsPack...>> : std::integral_constant<std::size_t, sizeof...(ExtentsPack)> {};

template <typename T>
static constexpr std::size_t extents_rank_v = extents_rank<std::remove_reference_t<T>>::value;

template <typename BaseExtents, std::size_t Dim, std::size_t N, std::size_t... I>
auto static_extent_transform(std::index_sequence<I...>)
    -> extents<typename BaseExtents::index_type, (I == Dim ? N : BaseExtents::static_extent(I))...>;

template <typename BaseExtents, typename Spec>
struct static_extents;

template <typename BaseExtents, std::size_t Dim, std::size_t N>
struct static_extents<BaseExtents, extent_with_static_dim<Dim, N>> {
    static_assert(Dim < BaseExtents::rank(), "extent_with_static_dim dimension must be smaller than extents rank.");
    using type = decltype(static_extent_transform<BaseExtents, Dim, N>(std::make_index_sequence<BaseExtents::rank()>{}));
};

template <typename BaseExtents, std::size_t N>
struct static_extents<BaseExtents, extent_with_static_last_dim<N>> {
    static_assert(BaseExtents::rank() > 0, "extent_with_static_last_dim requires rank greater than zero.");
    using type = typename static_extents<BaseExtents, extent_with_static_dim<BaseExtents::rank() - 1, N>>::type;
};

template <typename BaseExtents, typename Spec>
using static_extents_t = typename static_extents<BaseExtents, Spec>::type;

template <typename View, typename Extents, std::size_t... I>
Extents make_mdspan_extents(const View& view, std::index_sequence<I...>) {
    return Extents{static_cast<typename Extents::index_type>(mdspan_introspection_detail::extent(view, I))...};
}

template <typename View, std::size_t Rank, std::size_t... I>
std::array<std::size_t, Rank> make_mdspan_strides(const View& view, std::index_sequence<I...>) {
    return std::array<std::size_t, Rank>{static_cast<std::size_t>(mdspan_introspection_detail::stride(view, I))...};
}

template <typename Extents, std::size_t... I>
std::array<std::size_t, Extents::rank()> make_contiguous_strides(const Extents& extents, std::index_sequence<I...>) {
    constexpr std::size_t Rank = Extents::rank();
    std::array<std::size_t, Rank> strides{};
    std::size_t stride{1};
    for (std::size_t dimension = Rank; dimension > 0; --dimension) {
        strides[dimension - 1] = stride;
        stride *= static_cast<std::size_t>(extents.extent(dimension - 1));
    }
    return std::array<std::size_t, Rank>{strides[I]...};
}

template <typename ElementType>
inline constexpr std::size_t default_accessor_alignment_v = alignof(ElementType);

template <typename T, typename = std::void_t<> >
struct is_accessor_policy : std::false_type {};

template <typename T>
struct is_accessor_policy<T, std::void_t<typename T::template type<int>>> : std::true_type {};

template <typename T>
static constexpr bool is_accessor_policy_v = is_accessor_policy<T>::value;

template <typename AccessorPolicy, typename Value>
using accessor_policy_t = typename AccessorPolicy::template type<Value>;

template <typename T, typename = std::void_t<> >
struct is_accessor_type : std::false_type {};

template <typename T>
struct is_accessor_type<T, std::void_t<
    typename T::element_type,
    typename T::reference,
    typename T::data_handle_type>> : std::true_type {};

template <typename T>
static constexpr bool is_accessor_type_v = is_accessor_type<T>::value;

template <typename T>
static constexpr bool is_accessor_or_policy_v = is_accessor_policy_v<T> || is_accessor_type_v<T>;

template <typename AccessorOrPolicy, typename Value, bool = is_accessor_policy_v<AccessorOrPolicy>>
struct accessor_or_policy {
    using type = AccessorOrPolicy;
};

template <typename AccessorPolicy, typename Value>
struct accessor_or_policy<AccessorPolicy, Value, true> {
    using type = accessor_policy_t<AccessorPolicy, Value>;
};

template <typename AccessorOrPolicy, typename Value>
using accessor_or_policy_t = typename accessor_or_policy<AccessorOrPolicy, Value>::type;

template<typename Extents, typename Layout, typename Accessor, typename DataHandle, typename StridesFactory>
auto make_mdspan_from_extents(DataHandle data, Extents extents, StridesFactory strides_factory) {
    if constexpr (std::is_same_v<Layout, layout_right>) {
        typename layout_right::template mapping<Extents> mapping{extents};
        using Mdspan = mdspan<typename Accessor::element_type, Extents, layout_right, Accessor>;
        return Mdspan{data, mapping, Accessor{}};
    }
    else if constexpr (std::is_same_v<Layout, layout_stride>) {
        typename layout_stride::template mapping<Extents> mapping{extents, strides_factory()};
        using Mdspan = mdspan<typename Accessor::element_type, Extents, layout_stride, Accessor>;
        return Mdspan{data, mapping, Accessor{}};
    }
    else {
        static_assert(mdspan_introspection_detail::always_false_v<Layout>, "make_mdspan() is only implemented for layout_right and layout_stride");
        return mdspan<typename Accessor::element_type, Extents, layout_right, Accessor>{data};
    }
}

template<typename Extents, typename Layout, typename Accessor, typename View>
auto make_mdspan_impl(View& view) {
    constexpr std::size_t Rank = view_rank_v<View>;
    using ActualExtents = extract_extents_t<Extents, Rank>;

    auto extents = make_mdspan_extents<View, ActualExtents>(view, std::make_index_sequence<Rank>{});
    return make_mdspan_from_extents<ActualExtents, Layout, Accessor>(
        mdspan_introspection_detail::data_handle(view), extents,
        [&view]() { return make_mdspan_strides<View, Rank>(view, std::make_index_sequence<Rank>{}); });
}

template<typename Layout, typename Accessor, typename InputExtents, typename View>
auto make_mdspan_impl(View& view, InputExtents input_shape) {
    constexpr std::size_t Rank = view_rank_v<View>;
    using Extents = extract_extents_t<InputExtents, Rank>;

    static_assert(Extents::rank() == Rank, "The passed shape rank must match the view rank.");

    Extents extents{input_shape};
    return make_mdspan_from_extents<Extents, Layout, Accessor>(
        mdspan_introspection_detail::data_handle(view), extents,
        [&view]() { return make_mdspan_strides<View, Rank>(view, std::make_index_sequence<Rank>{}); });
}

template<typename Extents, typename Layout, typename Accessor, typename Value>
auto make_mdspan_pointer_impl(Value* data, Extents extents) {
    return make_mdspan_from_extents<Extents, Layout, Accessor>(
        data, extents,
        [&extents]() { return make_contiguous_strides(extents, std::make_index_sequence<Extents::rank()>{}); });
}

template<typename Extents, typename Layout, typename Accessor, typename Container>
auto make_mdspan_container_impl(Container& container, Extents extents) {
    return make_mdspan_pointer_impl<Extents, Layout, Accessor>(container.data(), extents);
}

}  // namespace make_mdspan_helpers
}  // namespace atlas
