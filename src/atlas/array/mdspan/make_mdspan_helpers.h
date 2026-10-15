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

#include <array>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "atlas/array/mdspan/mdspan.h"
#include "atlas/array/introspection.h"

namespace atlas::array {
template <typename ElementType, int Rank> class ArrayView;
template <typename ElementType, int Rank> class IndexView;
template <typename ElementType, int Rank> class LocalView;
}

namespace atlas::array::make_mdspan_helpers {

template <typename View>
struct is_non_owning_view : std::false_type {};

template <typename Element, int Rank>
struct is_non_owning_view<ArrayView<Element, Rank>> : std::true_type {};

template <typename Element, int Rank>
struct is_non_owning_view<IndexView<Element, Rank>> : std::true_type {};

template <typename Element, int Rank>
struct is_non_owning_view<LocalView<Element, Rank>> : std::true_type {};

template <typename Element, typename Extents, typename Layout, typename Accessor>
struct is_non_owning_view<mdspan<Element, Extents, Layout, Accessor>> : std::true_type {};

template <typename View>
using enable_rvalue_view_t = std::enable_if_t<
    !std::is_lvalue_reference_v<View> &&
    is_non_owning_view<std::remove_cv_t<std::remove_reference_t<View>>>::value, int>;

template <typename View>
using element_t = typename introspection::element_t<View>;

template <typename View>
using extents_t = typename introspection::extents_t<View>;

template <typename View>
using layout_t = typename introspection::layout_t<View>;

template <typename View>
using accessor_t = typename introspection::accessor_t<View>;

template <typename Container, typename = std::void_t<> >
struct is_container_like : std::false_type {};

template <typename Container>
struct is_container_like<Container, std::void_t<
    decltype(std::declval<Container&>().data()),
    decltype(std::declval<Container&>().size())>>
    : std::integral_constant<bool, !introspection::has_rank<Container>()> {};

template <typename Container>
static constexpr bool is_container_like_v = is_container_like<Container>::value;

template <typename Container>
struct container_descriptor {
    using extents_type = dextents<std::size_t, 1>;

    template <typename Source>
    static auto data(Source& source) { return source.data(); }

    static extents_type extents(const Container& source) { return extents_type{source.size()}; }
};

template <typename Element, std::size_t N, typename Allocator>
struct container_descriptor<std::vector<std::array<Element, N>, Allocator>> {
    using container_type = std::vector<std::array<Element, N>, Allocator>;
    using extents_type = atlas::array::extents<std::size_t, dynamic_extent, N>;

    static Element* data(container_type& source) { return reinterpret_cast<Element*>(source.data()); }
    static const Element* data(const container_type& source) { return reinterpret_cast<const Element*>(source.data()); }

    static extents_type extents(const container_type& source) { return extents_type{source.size()}; }
};

template <typename Container>
using container_descriptor_t = container_descriptor<std::remove_cv_t<Container>>;

template <typename Container>
using container_element_t = std::remove_pointer_t<decltype(container_descriptor_t<Container>::data(std::declval<Container&>()))>;

template <typename Container>
using enable_rvalue_container_t = std::enable_if_t<
    !std::is_lvalue_reference_v<Container> &&
    is_container_like_v<std::remove_reference_t<Container>>, int>;

template <typename T, std::size_t Rank>
struct to_extents { using type = T; };

template <typename T, std::size_t N, std::size_t Rank>
struct to_extents<std::array<T, N>, Rank> {
    using type = dextents<T, Rank>;
};

template <typename IndexType, std::size_t... ExtentsPack, std::size_t Rank>
struct to_extents<extents<IndexType, ExtentsPack...>, Rank> {
    using type = extents<IndexType, ExtentsPack...>;
};

template <typename T, std::size_t Rank>
using to_extents_t = typename to_extents<T, Rank>::type;

template <typename T, typename = std::void_t<>>
struct is_extent_like : std::false_type {};

template <typename T>
struct is_extent_like<T, std::void_t<decltype(std::declval<T>().rank())>> : std::true_type {};

template <typename T>
static constexpr bool is_extent_like_v = is_extent_like<T>::value;

template <typename T>
struct is_static_extent_spec : std::false_type {};

template <std::size_t Dim, std::size_t N>
struct is_static_extent_spec<fixed_dim<Dim, N>> : std::true_type {};

template <std::size_t N>
struct is_static_extent_spec<fixed_last_dim<N>> : std::true_type {};

template <typename T>
static constexpr bool is_static_extent_spec_v = is_static_extent_spec<T>::value;

namespace detail {
template <typename T, typename = std::void_t<> >
struct extents_rank;

template <typename T, std::size_t N>
struct extents_rank<std::array<T, N>> : std::integral_constant<std::size_t, N> {};

template <typename IndexType, std::size_t... ExtentsPack>
struct extents_rank<extents<IndexType, ExtentsPack...>> : std::integral_constant<std::size_t, sizeof...(ExtentsPack)> {};

}

template <typename Extents>
static constexpr std::size_t rank_from_extents() { return detail::extents_rank<std::remove_reference_t<Extents>>::value; }

template <typename BaseExtents, std::size_t Dim, std::size_t N, std::size_t... I>
auto static_extent_transform(std::index_sequence<I...>)
    -> extents<typename BaseExtents::index_type, (I == Dim ? N : BaseExtents::static_extent(I))...>;

template <typename BaseExtents, typename Spec>
struct static_extents;

template <typename BaseExtents, std::size_t Dim, std::size_t N>
struct static_extents<BaseExtents, fixed_dim<Dim, N>> {
    static_assert(Dim < BaseExtents::rank(), "fixed_dim dimension must be smaller than extents rank.");
    using type = decltype(static_extent_transform<BaseExtents, Dim, N>(std::make_index_sequence<BaseExtents::rank()>{}));
};

template <typename BaseExtents, std::size_t N>
struct static_extents<BaseExtents, fixed_last_dim<N>> {
    static_assert(BaseExtents::rank() > 0, "fixed_last_dim requires rank greater than zero.");
    using type = typename static_extents<BaseExtents, fixed_dim<BaseExtents::rank() - 1, N>>::type;
};

template <typename BaseExtents, typename Spec>
using static_extents_t = typename static_extents<BaseExtents, Spec>::type;

template <typename View, typename Extents, std::size_t... I>
Extents make_mdspan_extents(const View& view, std::index_sequence<I...>) {
    return Extents{static_cast<typename Extents::index_type>(extent<I>(view))...};
}

template <typename View, std::size_t Rank, std::size_t... I>
std::array<std::size_t, Rank> make_mdspan_strides(const View& view, std::index_sequence<I...>) {
    [[maybe_unused]] auto checked_stride = [](std::size_t dimension, auto value) {
        if constexpr (std::is_signed_v<decltype(value)>) {
            if (value < 0) {
                throw std::invalid_argument("make_mdspan(): negative source stride at dimension " +
                                            std::to_string(dimension) + ": " + std::to_string(value));
            }
        }
        return static_cast<std::size_t>(value);
    };
    return std::array<std::size_t, Rank>{checked_stride(I, stride<I>(view))...};
}

template <typename Extents, std::size_t... I>
std::array<std::size_t, Extents::rank()> make_contiguous_strides(const Extents& extents, std::index_sequence<I...>) {
    constexpr std::size_t Rank = Extents::rank();
    std::array<std::size_t, Rank> strides{};
    std::size_t stride{1};
    for (std::size_t dimension = Rank; dimension > 0; --dimension) {
        strides[dimension - 1] = stride;
        const auto extent = static_cast<std::size_t>(extents.extent(dimension - 1));
        stride *= extent == 0 ? 1 : extent;
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

// A type is treated as a layout policy when it is none of the other recognized make_mdspan
// template categories (extents, static extent spec, accessor policy, or accessor type). This
// covers the mdspan layout policies such as layout_right, layout_stride, and layout_left.
template <typename T>
static constexpr bool is_layout_like_v = !is_extent_like_v<T> && !is_static_extent_spec_v<T> &&
                                         !is_accessor_policy_v<T> && !is_accessor_type_v<T>;

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

template <typename ValueAt>
std::string format_dimensions(std::size_t rank, ValueAt value_at) {
    std::ostringstream result;
    result << "{";
    for (std::size_t dimension = 0; dimension < rank; ++dimension) {
        result << (dimension == 0 ? "" : ", ") << value_at(dimension);
    }
    result << "}";
    return result.str();
}

template<typename Extents, typename Layout, typename Accessor, typename DataHandle, typename StridesFactory>
auto make_mdspan_from_extents(DataHandle data, Extents extents, StridesFactory strides_factory,
                             Accessor accessor = Accessor{}) {
    if constexpr (std::is_same_v<Layout, layout_right> || std::is_same_v<Layout, layout_left>) {
        typename Layout::template mapping<Extents> mapping{extents};
        const auto strides = strides_factory();
        if (!introspection::valid_contiguous_layout(mapping, strides)) {
            const auto layout_name = std::is_same_v<Layout, layout_left> ? "layout_left" : "layout_right";
            std::ostringstream message;
            message << "make_mdspan(): source strides are incompatible with " << layout_name << "; extents="
                << format_dimensions(Extents::rank(), [&](auto index) { return extents.extent(index); })
                << ", source strides="
                << format_dimensions(Extents::rank(), [&](auto index) { return strides[index]; })
                << ", required strides="
                << format_dimensions(Extents::rank(), [&](auto index) { return mapping.stride(index); });
            throw std::invalid_argument(message.str());
        }
        using Mdspan = mdspan<typename Accessor::element_type, Extents, Layout, Accessor>;
        return Mdspan{data, mapping, accessor};
    }
    else if constexpr (std::is_same_v<Layout, layout_stride>) {
        const auto strides = strides_factory();
        std::array<std::size_t, Extents::rank()> shape{};
        for (std::size_t dimension = 0; dimension < Extents::rank(); ++dimension) {
            shape[dimension] = extents.extent(dimension);
        }
        if (!introspection::valid_layout_stride(shape, strides)) {
            std::ostringstream message;
            message << "make_mdspan(): layout_stride requires positive, non-overlapping strides; extents="
                    << format_dimensions(Extents::rank(), [&](auto dimension) { return shape[dimension]; })
                    << ", strides="
                    << format_dimensions(Extents::rank(), [&](auto dimension) { return strides[dimension]; });
            throw std::invalid_argument(message.str());
        }
        typename layout_stride::template mapping<Extents> mapping{extents, strides};
        using Mdspan = mdspan<typename Accessor::element_type, Extents, layout_stride, Accessor>;
        return Mdspan{data, mapping, accessor};
    }
    else {
        static_assert(introspection::always_false_v<Layout>, "make_mdspan() is only implemented for layout_left, layout_right and layout_stride");
        return mdspan<typename Accessor::element_type, Extents, layout_right, Accessor>{data};
    }
}

template <typename View, typename = void>
struct has_view_accessor : std::false_type {};

template <typename View>
struct has_view_accessor<View, std::void_t<typename View::accessor_type,
                                         decltype(std::declval<const View&>().accessor())>> : std::true_type {};

template <typename Accessor, typename View>
auto make_view_accessor(const View& view) {
    if constexpr (has_view_accessor<View>::value) {
        if constexpr (std::is_same_v<Accessor, typename View::accessor_type>) {
            return view.accessor();
        }
        else {
            return Accessor{};
        }
    }
    else {
        return Accessor{};
    }
}

template<typename Extents, typename Layout, typename Accessor, typename View>
auto make_mdspan_impl(View& view) {
    using SourceLayout = layout_t<View>;
    static_assert((!(std::is_same_v<SourceLayout, layout_right> && std::is_same_v<Layout, layout_left>) &&
                  !(std::is_same_v<SourceLayout, layout_left> && std::is_same_v<Layout, layout_right>)),
                 "make_mdspan(): cannot convert between declared left and right layouts");
    constexpr std::size_t Rank = rank<View>();
    using ActualExtents = to_extents_t<Extents, Rank>;

    static_assert(ActualExtents::rank() == Rank, "The requested extents rank must match the view rank.");

    auto extents = make_mdspan_extents<View, ActualExtents>(view, std::make_index_sequence<Rank>{});
    return make_mdspan_from_extents<ActualExtents, Layout, Accessor>(
        data_handle(view), extents,
        [&view]() { return make_mdspan_strides<View, Rank>(view, std::make_index_sequence<Rank>{}); },
        make_view_accessor<Accessor>(view));
}

template<typename Layout, typename Accessor, typename InputExtents, typename View>
auto make_mdspan_impl(View& view, InputExtents input_shape) {
    using SourceLayout = layout_t<View>;
    static_assert((!(std::is_same_v<SourceLayout, layout_right> && std::is_same_v<Layout, layout_left>) &&
                  !(std::is_same_v<SourceLayout, layout_left> && std::is_same_v<Layout, layout_right>)),
                 "make_mdspan(): cannot convert between declared left and right layouts");
    constexpr std::size_t Rank = rank<View>();
    using Extents = to_extents_t<InputExtents, Rank>;

    static_assert(Extents::rank() == Rank, "The passed shape rank must match the view rank.");

    Extents extents{input_shape};
    return make_mdspan_from_extents<Extents, Layout, Accessor>(
        data_handle(view), extents,
        [&view]() { return make_mdspan_strides<View, Rank>(view, std::make_index_sequence<Rank>{}); },
        make_view_accessor<Accessor>(view));
}

template<typename Extents, typename Layout, typename Accessor, typename Value>
auto make_mdspan_pointer_impl(Value* data, Extents extents) {
    return make_mdspan_from_extents<Extents, Layout, Accessor>(
        data, extents,
        [&extents]() {
            std::array<std::size_t, Extents::rank()> strides{};
            if constexpr (std::is_same_v<Layout, layout_left>) {
                typename layout_left::template mapping<Extents> mapping{extents};
                for (std::size_t dimension = 0; dimension < Extents::rank(); ++dimension) {
                    strides[dimension] = mapping.stride(dimension);
                }
            }
            else {
                strides = make_contiguous_strides(extents, std::make_index_sequence<Extents::rank()>{});
            }
            return strides;
        });
}

template<typename Extents, typename Layout, typename Accessor, typename Container>
auto make_mdspan_container_impl(Container& container, Extents extents) {
    return make_mdspan_pointer_impl<Extents, Layout, Accessor>(container_descriptor_t<Container>::data(container), extents);
}

}  // namespace atlas::array::make_mdspan_helpers
