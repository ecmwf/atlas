#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "atlas/mdspan.h"

namespace atlas {
namespace array {

template <std::size_t Dim, std::size_t N>
struct StaticExtent {};

template <std::size_t N>
struct StaticLastExtent {};

namespace detail {

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

template <typename Container, typename = std::void_t<> >
struct is_container_like : std::false_type {};

template <typename Container>
struct is_container_like<Container, std::void_t<
    decltype(std::declval<Container&>().data()),
    decltype(std::declval<Container&>().size())>>
    : std::integral_constant<bool, !has_view_rank_v<Container>> {};

template <typename Container>
static constexpr bool is_container_like_v = is_container_like<Container>::value;

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

template <typename View, typename = std::void_t<>>
struct mdspan_value {
    using type = std::remove_pointer_t<std::remove_reference_t<decltype(std::declval<View&>().data())>>;
};

template <typename View>
struct mdspan_value<View, std::void_t<decltype(std::declval<View&>().data_handle())>> {
    using type = std::remove_pointer_t<std::remove_reference_t<decltype(std::declval<View&>().data_handle())>>;
};

template <typename View>
using mdspan_value_t = typename mdspan_value<View>::type;

template <typename View, typename = std::void_t<> >
struct view_accessor {
    using type = restrict_accessor<mdspan_value_t<View>>;
};

template <typename View>
struct view_accessor<View, std::void_t<typename view_type_t<View>::accessor_type>> {
    using type = typename view_type_t<View>::accessor_type;
};

template <typename View>
using view_accessor_t = typename view_accessor<View>::type;

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

template <typename View, typename Extents, std::size_t... I>
Extents make_mdspan_extents(const View& view, std::index_sequence<I...>) {
    return Extents{static_cast<typename Extents::index_type>(extent(view, I))...};
}

template <typename View, std::size_t Rank, std::size_t... I>
std::array<std::size_t, Rank> make_mdspan_strides(const View& view, std::index_sequence<I...>) {
    return std::array<std::size_t, Rank>{static_cast<std::size_t>(stride(view, I))...};
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

template <typename>
[[maybe_unused]] inline static constexpr bool always_false_v = false;

template <typename ElementType>
inline constexpr std::size_t default_accessor_alignment_v = alignof(ElementType);

template <template <typename, std::size_t> typename AccessorPolicy>
struct bind_default_alignment {
    template <typename ElementType>
    using type = AccessorPolicy<ElementType, default_accessor_alignment_v<ElementType>>;
};

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
    using Value = std::remove_cv_t<mdspan_value_t<View>>;
    return (((static_cast<std::size_t>(stride(view, I)) * sizeof(Value)) % alignment == 0) && ...);
}

template <std::size_t Dim, typename View, std::size_t... I>
bool dimension_start_strides_are_aligned(const View& view, std::size_t alignment, std::index_sequence<I...>) {
    return dimension_start_strides_are_aligned(
        view, alignment, std::index_sequence<(I < Dim ? I : I + 1)...>{});
}

template<typename Extents, typename Layout, typename Accessor, typename View>
auto make_mdspan_impl(View& view) {
    using Value = mdspan_value_t<View>;
    constexpr std::size_t Rank = view_rank_v<View>;
    using ActualExtents = extract_extents_t<Extents, Rank>;

    auto extents = make_mdspan_extents<View, ActualExtents>(view, std::make_index_sequence<Rank>{});
    if constexpr (std::is_same_v<Layout, layout_right>) {
        typename layout_right::template mapping<ActualExtents> mapping{extents};
        return mdspan<Value, ActualExtents, layout_right, Accessor>(data_handle(view), mapping, Accessor{});
    }
    else if constexpr (std::is_same_v<Layout, layout_stride>) {
        typename layout_stride::template mapping<ActualExtents> mapping{
            extents, make_mdspan_strides<View, Rank>(view, std::make_index_sequence<Rank>{})};
        return mdspan<Value, ActualExtents, layout_stride, Accessor>(data_handle(view), mapping, Accessor{});
    }
    else {
        static_assert(always_false_v<Layout>, "make_mdspan() is only implemented for layout_right and layout_stride");
    }
}

template<typename Layout, typename Accessor, typename InputExtents, typename View>
auto make_mdspan_impl(View& view, InputExtents input_shape) {
    using Value = mdspan_value_t<View>;
    constexpr std::size_t Rank = view_rank_v<View>;
    using Extents = extract_extents_t<InputExtents, Rank>;

    static_assert(Extents::rank() == Rank, "The passed shape rank must match the view rank.");

    Extents extents{input_shape};
    if constexpr (std::is_same_v<Layout, layout_right>) {
        typename layout_right::template mapping<Extents> mapping{extents};
        return mdspan<Value, Extents, layout_right, Accessor>(data_handle(view), mapping, Accessor{});
    }
    else if constexpr (std::is_same_v<Layout, layout_stride>) {
        typename layout_stride::template mapping<Extents> mapping{
            extents, make_mdspan_strides<View, Rank>(view, std::make_index_sequence<Rank>{})};
        return mdspan<Value, Extents, layout_stride, Accessor>(data_handle(view), mapping, Accessor{});
    }
    else {
        static_assert(always_false_v<Layout>, "make_mdspan() is only implemented for layout_right and layout_stride");
    }
}

template<typename Extents, typename Layout, typename Accessor, typename Container>
auto make_mdspan_container_impl(Container& container, Extents extents) {
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;

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
        static_assert(always_false_v<Layout>, "make_mdspan() is only implemented for layout_right and layout_stride");
    }
}

template<typename Extents, typename Layout, typename Accessor, typename Value>
auto make_mdspan_pointer_impl(Value* data, Extents extents) {
    if constexpr (std::is_same_v<Layout, layout_right>) {
        typename layout_right::template mapping<Extents> mapping{extents};
        return mdspan<Value, Extents, layout_right, Accessor>(data, mapping, Accessor{});
    }
    else if constexpr (std::is_same_v<Layout, layout_stride>) {
        typename layout_stride::template mapping<Extents> mapping{
            extents, make_contiguous_strides(extents, std::make_index_sequence<Extents::rank()>{})};
        return mdspan<Value, Extents, layout_stride, Accessor>(data, mapping, Accessor{});
    }
    else {
        static_assert(always_false_v<Layout>, "make_mdspan() is only implemented for layout_right and layout_stride");
    }
}

}  // namespace detail

template <typename View>
bool can_use_layout_right(const View& view) {
    constexpr std::size_t Rank = detail::view_rank_v<View>;
    if constexpr (Rank <= 1) {
        return true;
    }
    else {
        return detail::has_layout_right_strides(view, std::make_index_sequence<Rank>{});
    }
}

template <typename View>
bool is_aligned(View& view, std::size_t alignment) {
    if (alignment == 0) {
        return false;
    }
    auto address = reinterpret_cast<std::uintptr_t>(detail::data_handle(view));
    return address % alignment == 0;
}

template <std::size_t Dim, typename View>
bool is_dimension_aligned(View& view, std::size_t alignment) {
    constexpr std::size_t Rank = detail::view_rank_v<View>;

    static_assert(Dim < Rank, "Dimension must be smaller than view rank.");

    if (!is_aligned(view, alignment)) {
        return false;
    }

    if constexpr (Rank <= 1) {
        return true;
    }
    else {
        return detail::dimension_start_strides_are_aligned<Dim>(view, alignment, std::make_index_sequence<Rank - 1>{});
    }
}

template <typename View>
bool is_last_dimension_aligned(View& view, std::size_t alignment) {
    constexpr std::size_t Rank = detail::view_rank_v<View>;
    return is_dimension_aligned<Rank - 1>(view, alignment);
}

template<
    typename Extents,
    typename Layout = layout_stride,
    template <typename> typename AccessorPolicy = restrict_accessor,
    typename View,
    typename std::enable_if_t<detail::is_extent_like_v<Extents>, int> = 0
>
auto make_mdspan(View& view) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    return detail::make_mdspan_impl<Extents, Layout, Accessor>(view);
}

template<
    typename Layout,
    typename View,
    typename std::enable_if_t<!detail::is_extent_like_v<Layout> && !detail::is_static_extent_spec_v<Layout> && !detail::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    return detail::make_mdspan_impl<detail::view_extents_t<View>, Layout, detail::view_accessor_t<View>>(view);
}

template<
    typename Layout,
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<!detail::is_extent_like_v<Layout> && !detail::is_static_extent_spec_v<Layout> && !detail::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    return detail::make_mdspan_impl<detail::view_extents_t<View>, Layout, Accessor>(view);
}

template<
    typename StaticExtentsSpec,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using Extents = detail::static_extents_t<detail::view_extents_t<View>, StaticExtentsSpec>;
    return detail::make_mdspan_impl<Extents, detail::view_layout_t<View>, detail::view_accessor_t<View>>(view);
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec> && !detail::is_extent_like_v<Layout>, int> = 0
>
auto make_mdspan(View& view) {
    using Extents = detail::static_extents_t<detail::view_extents_t<View>, StaticExtentsSpec>;
    return detail::make_mdspan_impl<Extents, Layout, detail::view_accessor_t<View>>(view);
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    using Extents = detail::static_extents_t<detail::view_extents_t<View>, StaticExtentsSpec>;
    return detail::make_mdspan_impl<Extents, Layout, Accessor>(view);
}

template<
    typename StaticExtentsSpec,
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    using Extents = detail::static_extents_t<detail::view_extents_t<View>, StaticExtentsSpec>;
    return detail::make_mdspan_impl<Extents, detail::view_layout_t<View>, Accessor>(view);
}

template<
    typename StaticExtentsSpec,
    template <typename, std::size_t> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = typename detail::bind_default_alignment<AccessorPolicy>::template type<Value>;
    using Extents = detail::static_extents_t<detail::view_extents_t<View>, StaticExtentsSpec>;
    return detail::make_mdspan_impl<Extents, detail::view_layout_t<View>, Accessor>(view);
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    template <typename, std::size_t> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = typename detail::bind_default_alignment<AccessorPolicy>::template type<Value>;
    using Extents = detail::static_extents_t<detail::view_extents_t<View>, StaticExtentsSpec>;
    return detail::make_mdspan_impl<Extents, Layout, Accessor>(view);
}

template<
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<!detail::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    return make_mdspan<detail::view_extents_t<View>, detail::view_layout_t<View>, AccessorPolicy>(view);
}

template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<!detail::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    return make_mdspan<detail::view_extents_t<View>, detail::view_layout_t<View>,
                       detail::bind_default_alignment<AccessorPolicy>::template type>(view);
}

template<typename View, typename std::enable_if_t<!detail::is_container_like_v<View>, int> = 0>
auto make_mdspan(View& view) {
    return detail::make_mdspan_impl<detail::view_extents_t<View>, detail::view_layout_t<View>,
                                    detail::view_accessor_t<View>>(view);
}

template<
    typename Layout = layout_right,
    template <typename> typename AccessorPolicy = restrict_accessor,
    typename Container,
    typename std::enable_if_t<detail::is_container_like_v<Container> && !detail::is_extent_like_v<Layout>, int> = 0
>
auto make_mdspan(Container& container) {
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    using Extents = dextents<std::size_t, 1>;
    using Accessor = AccessorPolicy<Value>;
    return detail::make_mdspan_container_impl<Extents, Layout, Accessor>(container, Extents{container.size()});
}

template<
    template <typename> typename AccessorPolicy,
    typename Container,
    typename std::enable_if_t<detail::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    return make_mdspan<layout_right, AccessorPolicy>(container);
}

template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename Container,
    typename std::enable_if_t<detail::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    return make_mdspan<layout_right, detail::bind_default_alignment<AccessorPolicy>::template type>(container);
}

template<
    typename Layout,
    template <typename> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!detail::is_extent_like_v<Layout> && !detail::is_static_extent_spec_v<Layout> && !detail::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    return detail::make_mdspan_impl<Layout, Accessor>(view, input_shape);
}

template<
    typename Layout,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!detail::is_extent_like_v<Layout> && !detail::is_static_extent_spec_v<Layout> && !detail::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    return detail::make_mdspan_impl<Layout, detail::view_accessor_t<View>>(view, input_shape);
}

template<
    typename StaticExtentsSpec,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    constexpr std::size_t Rank = detail::view_rank_v<View>;
    using BaseExtents = detail::extract_extents_t<InputExtents, Rank>;
    using Extents = detail::static_extents_t<BaseExtents, StaticExtentsSpec>;
    return detail::make_mdspan_impl<detail::view_layout_t<View>, detail::view_accessor_t<View>>(view, Extents{input_shape});
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec> && !detail::is_extent_like_v<Layout>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    constexpr std::size_t Rank = detail::view_rank_v<View>;
    using BaseExtents = detail::extract_extents_t<InputExtents, Rank>;
    using Extents = detail::static_extents_t<BaseExtents, StaticExtentsSpec>;
    return detail::make_mdspan_impl<Layout, detail::view_accessor_t<View>>(view, Extents{input_shape});
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    template <typename> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    constexpr std::size_t Rank = detail::view_rank_v<View>;
    using BaseExtents = detail::extract_extents_t<InputExtents, Rank>;
    using Extents = detail::static_extents_t<BaseExtents, StaticExtentsSpec>;
    return detail::make_mdspan_impl<Layout, Accessor>(view, Extents{input_shape});
}

template<
    typename StaticExtentsSpec,
    template <typename> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    constexpr std::size_t Rank = detail::view_rank_v<View>;
    using BaseExtents = detail::extract_extents_t<InputExtents, Rank>;
    using Extents = detail::static_extents_t<BaseExtents, StaticExtentsSpec>;
    return detail::make_mdspan_impl<detail::view_layout_t<View>, Accessor>(view, Extents{input_shape});
}

template<
    typename StaticExtentsSpec,
    template <typename, std::size_t> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = typename detail::bind_default_alignment<AccessorPolicy>::template type<Value>;
    constexpr std::size_t Rank = detail::view_rank_v<View>;
    using BaseExtents = detail::extract_extents_t<InputExtents, Rank>;
    using Extents = detail::static_extents_t<BaseExtents, StaticExtentsSpec>;
    return detail::make_mdspan_impl<detail::view_layout_t<View>, Accessor>(view, Extents{input_shape});
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    template <typename, std::size_t> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<detail::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using Value = detail::mdspan_value_t<View>;
    using Accessor = typename detail::bind_default_alignment<AccessorPolicy>::template type<Value>;
    constexpr std::size_t Rank = detail::view_rank_v<View>;
    using BaseExtents = detail::extract_extents_t<InputExtents, Rank>;
    using Extents = detail::static_extents_t<BaseExtents, StaticExtentsSpec>;
    return detail::make_mdspan_impl<Layout, Accessor>(view, Extents{input_shape});
}

template<
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!detail::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    return detail::make_mdspan_impl<detail::view_layout_t<View>, detail::view_accessor_t<View>>(view, input_shape);
}

template<
    template <typename> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!detail::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    return make_mdspan<detail::view_layout_t<View>, AccessorPolicy>(view, input_shape);
}

template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!detail::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    return make_mdspan<detail::view_layout_t<View>, detail::bind_default_alignment<AccessorPolicy>::template type>(view, input_shape);
}

template<
    typename Layout = layout_right,
    template <typename> typename AccessorPolicy = restrict_accessor,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<detail::is_container_like_v<Container> && !detail::is_extent_like_v<Layout>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    constexpr std::size_t Rank = detail::extents_rank_v<InputExtents>;
    using Extents = detail::extract_extents_t<InputExtents, Rank>;
    using Accessor = AccessorPolicy<Value>;
    return detail::make_mdspan_container_impl<Extents, Layout, Accessor>(container, Extents{input_shape});
}

template<
    template <typename> typename AccessorPolicy,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<detail::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    return make_mdspan<layout_right, AccessorPolicy>(container, input_shape);
}

template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<detail::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    return make_mdspan<layout_right, detail::bind_default_alignment<AccessorPolicy>::template type>(container, input_shape);
}

template<
    typename Layout = layout_right,
    template <typename> typename AccessorPolicy = restrict_accessor,
    typename Value,
    typename InputExtents,
    typename std::enable_if_t<!detail::is_extent_like_v<Layout>, int> = 0
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    constexpr std::size_t Rank = detail::extents_rank_v<InputExtents>;
    using Extents = detail::extract_extents_t<InputExtents, Rank>;
    using Accessor = AccessorPolicy<Value>;
    return detail::make_mdspan_pointer_impl<Extents, Layout, Accessor>(data, Extents{input_shape});
}

template<
    template <typename> typename AccessorPolicy,
    typename Value,
    typename InputExtents
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    return make_mdspan<layout_right, AccessorPolicy>(data, input_shape);
}

template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename Value,
    typename InputExtents
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    return make_mdspan<layout_right, detail::bind_default_alignment<AccessorPolicy>::template type>(data, input_shape);
}

}  // namespace array
}  // namespace atlas
