#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "atlas/mdspan.h"
#include "atlas/mdspan_introspection.h"

namespace atlas {
namespace make_mdspan_helpers {

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
struct is_static_extent_spec<StaticExtent<Dim, N>> : std::true_type {};

template <std::size_t N>
struct is_static_extent_spec<StaticLastExtent<N>> : std::true_type {};

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
struct static_extents<BaseExtents, StaticExtent<Dim, N>> {
    static_assert(Dim < BaseExtents::rank(), "StaticExtent dimension must be smaller than extents rank.");
    using type = decltype(static_extent_transform<BaseExtents, Dim, N>(std::make_index_sequence<BaseExtents::rank()>{}));
};

template <typename BaseExtents, std::size_t N>
struct static_extents<BaseExtents, StaticLastExtent<N>> {
    static_assert(BaseExtents::rank() > 0, "StaticLastExtent requires rank greater than zero.");
    using type = typename static_extents<BaseExtents, StaticExtent<BaseExtents::rank() - 1, N>>::type;
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

template <template <typename, std::size_t> typename AccessorPolicy>
struct bind_default_alignment {
    template <typename ElementType>
    using type = AccessorPolicy<ElementType, default_accessor_alignment_v<ElementType>>;
};

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

template<typename Extents, typename Layout, typename Accessor, typename View>
auto make_mdspan_impl(View& view) {
    using Value = typename Accessor::element_type;
    constexpr std::size_t Rank = view_rank_v<View>;
    using ActualExtents = extract_extents_t<Extents, Rank>;

    auto extents = make_mdspan_extents<View, ActualExtents>(view, std::make_index_sequence<Rank>{});
    if constexpr (std::is_same_v<Layout, layout_right>) {
        typename layout_right::template mapping<ActualExtents> mapping{extents};
        return mdspan<Value, ActualExtents, layout_right, Accessor>(mdspan_introspection_detail::data_handle(view), mapping, Accessor{});
    }
    else if constexpr (std::is_same_v<Layout, layout_stride>) {
        typename layout_stride::template mapping<ActualExtents> mapping{
            extents, make_mdspan_strides<View, Rank>(view, std::make_index_sequence<Rank>{})};
        return mdspan<Value, ActualExtents, layout_stride, Accessor>(mdspan_introspection_detail::data_handle(view), mapping, Accessor{});
    }
    else {
        static_assert(mdspan_introspection_detail::always_false_v<Layout>, "make_mdspan() is only implemented for layout_right and layout_stride");
    }
}

template<typename Layout, typename Accessor, typename InputExtents, typename View>
auto make_mdspan_impl(View& view, InputExtents input_shape) {
    using Value = typename Accessor::element_type;
    constexpr std::size_t Rank = view_rank_v<View>;
    using Extents = extract_extents_t<InputExtents, Rank>;

    static_assert(Extents::rank() == Rank, "The passed shape rank must match the view rank.");

    Extents extents{input_shape};
    if constexpr (std::is_same_v<Layout, layout_right>) {
        typename layout_right::template mapping<Extents> mapping{extents};
        return mdspan<Value, Extents, layout_right, Accessor>(mdspan_introspection_detail::data_handle(view), mapping, Accessor{});
    }
    else if constexpr (std::is_same_v<Layout, layout_stride>) {
        typename layout_stride::template mapping<Extents> mapping{
            extents, make_mdspan_strides<View, Rank>(view, std::make_index_sequence<Rank>{})};
        return mdspan<Value, Extents, layout_stride, Accessor>(mdspan_introspection_detail::data_handle(view), mapping, Accessor{});
    }
    else {
        static_assert(mdspan_introspection_detail::always_false_v<Layout>, "make_mdspan() is only implemented for layout_right and layout_stride");
    }
}

template<typename Extents, typename Layout, typename Accessor, typename Container>
auto make_mdspan_container_impl(Container& container, Extents extents) {
    using Value = typename Accessor::element_type;

    if constexpr (std::is_same_v<Layout, layout_right>) {
        typename layout_right::template mapping<Extents> mapping{extents};
        return mdspan<Value, Extents, layout_right, Accessor>(container.data(), mapping, Accessor{});
    }
    else if constexpr (std::is_same_v<Layout, layout_stride>) {
        typename layout_stride::template mapping<Extents> mapping{
            extents, make_contiguous_strides(extents, std::make_index_sequence<Extents::rank()>{})};
        return mdspan<Value, Extents, layout_stride, Accessor>(container.data(), mapping, Accessor{});
    }
    else {
        static_assert(mdspan_introspection_detail::always_false_v<Layout>, "make_mdspan() is only implemented for layout_right and layout_stride");
    }
}

template<typename Extents, typename Layout, typename Accessor, typename Value>
auto make_mdspan_pointer_impl(Value* data, Extents extents) {
    using ElementType = typename Accessor::element_type;
    if constexpr (std::is_same_v<Layout, layout_right>) {
        typename layout_right::template mapping<Extents> mapping{extents};
        return mdspan<ElementType, Extents, layout_right, Accessor>(data, mapping, Accessor{});
    }
    else if constexpr (std::is_same_v<Layout, layout_stride>) {
        typename layout_stride::template mapping<Extents> mapping{
            extents, make_contiguous_strides(extents, std::make_index_sequence<Extents::rank()>{})};
        return mdspan<ElementType, Extents, layout_stride, Accessor>(data, mapping, Accessor{});
    }
    else {
        static_assert(mdspan_introspection_detail::always_false_v<Layout>, "make_mdspan() is only implemented for layout_right and layout_stride");
    }
}

}  // namespace make_mdspan_helpers
}  // namespace atlas
