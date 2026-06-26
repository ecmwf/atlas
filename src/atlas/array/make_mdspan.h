#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "atlas/make_mdspan_helpers.h"
#include "atlas/mdspan.h"
#include "atlas/mdspan_introspection.h"

namespace atlas {

template<
    typename Extents,
    typename Layout = layout_stride,
    template <typename> typename AccessorPolicy = restrict_accessor,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_extent_like_v<Extents>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    return make_mdspan_impl<Extents, Layout, Accessor>(view);
}

template<
    typename Layout,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_static_extent_spec_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    return make_mdspan_impl<view_extents_t<View>, Layout, view_accessor_t<View>>(view);
}

template<
    typename Layout,
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_static_extent_spec_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    return make_mdspan_impl<view_extents_t<View>, Layout, Accessor>(view);
}

template<
    typename StaticExtentsSpec,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan_impl<Extents, view_layout_t<View>, view_accessor_t<View>>(view);
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec> && !make_mdspan_helpers::is_extent_like_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan_impl<Extents, Layout, view_accessor_t<View>>(view);
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan_impl<Extents, Layout, Accessor>(view);
}

template<
    typename StaticExtentsSpec,
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan_impl<Extents, view_layout_t<View>, Accessor>(view);
}

template<
    typename StaticExtentsSpec,
    template <typename, std::size_t> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = typename bind_default_alignment<AccessorPolicy>::template type<Value>;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan_impl<Extents, view_layout_t<View>, Accessor>(view);
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    template <typename, std::size_t> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = typename bind_default_alignment<AccessorPolicy>::template type<Value>;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan_impl<Extents, Layout, Accessor>(view);
}

template<
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    return make_mdspan<view_extents_t<View>, view_layout_t<View>, AccessorPolicy>(view);
}

template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    return make_mdspan<view_extents_t<View>, view_layout_t<View>,
                       bind_default_alignment<AccessorPolicy>::template type>(view);
}

template<
    typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_policy_v<AccessorPolicy> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = accessor_policy_t<AccessorPolicy, Value>;
    return make_mdspan_impl<view_extents_t<View>, view_layout_t<View>, Accessor>(view);
}

template<
    typename Accessor,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_type_v<Accessor> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    return make_mdspan_impl<view_extents_t<View>, view_layout_t<View>, Accessor>(view);
}

template<
    typename Layout,
    typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_static_extent_spec_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              make_mdspan_helpers::is_accessor_policy_v<AccessorPolicy> &&
                              !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = accessor_policy_t<AccessorPolicy, Value>;
    return make_mdspan_impl<view_extents_t<View>, Layout, Accessor>(view);
}

template<
    typename Layout,
    typename Accessor,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_static_extent_spec_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              make_mdspan_helpers::is_accessor_type_v<Accessor> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    return make_mdspan_impl<view_extents_t<View>, Layout, Accessor>(view);
}

template<typename View, typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    return make_mdspan_impl<view_extents_t<View>, view_layout_t<View>,
                                    view_accessor_t<View>>(view);
}

template<
    typename Layout = layout_right,
    template <typename> typename AccessorPolicy = restrict_accessor,
    typename Container,
    typename std::enable_if_t<make_mdspan_helpers::is_container_like_v<Container> && !make_mdspan_helpers::is_extent_like_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout>, int> = 0
>
auto make_mdspan(Container& container) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    using Extents = dextents<std::size_t, 1>;
    using Accessor = AccessorPolicy<Value>;
    return make_mdspan_container_impl<Extents, Layout, Accessor>(container, Extents{container.size()});
}

template<
    template <typename> typename AccessorPolicy,
    typename Container,
    typename std::enable_if_t<make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    return make_mdspan<layout_right, AccessorPolicy>(container);
}

template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename Container,
    typename std::enable_if_t<make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    using namespace make_mdspan_helpers;
    return make_mdspan<layout_right, bind_default_alignment<AccessorPolicy>::template type>(container);
}

template<
    typename AccessorPolicy,
    typename Container,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_policy_v<AccessorPolicy> && make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    using Extents = dextents<std::size_t, 1>;
    using Accessor = accessor_policy_t<AccessorPolicy, Value>;
    return make_mdspan_container_impl<Extents, layout_right, Accessor>(container, Extents{container.size()});
}

template<
    typename Layout,
    typename AccessorPolicy,
    typename Container,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_accessor_policy_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_type_v<Layout> && make_mdspan_helpers::is_accessor_policy_v<AccessorPolicy> &&
                              make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    using Extents = dextents<std::size_t, 1>;
    using Accessor = accessor_policy_t<AccessorPolicy, Value>;
    return make_mdspan_container_impl<Extents, Layout, Accessor>(container, Extents{container.size()});
}

template<
    typename Accessor,
    typename Container,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_type_v<Accessor> && make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    using namespace make_mdspan_helpers;
    using Extents = dextents<std::size_t, 1>;
    return make_mdspan_container_impl<Extents, layout_right, Accessor>(container, Extents{container.size()});
}

template<
    typename Layout,
    typename Accessor,
    typename Container,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_accessor_policy_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_type_v<Layout> && make_mdspan_helpers::is_accessor_type_v<Accessor> &&
                              make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    using namespace make_mdspan_helpers;
    using Extents = dextents<std::size_t, 1>;
    return make_mdspan_container_impl<Extents, Layout, Accessor>(container, Extents{container.size()});
}

template<
    typename Layout,
    template <typename> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_static_extent_spec_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    return make_mdspan_impl<Layout, Accessor>(view, input_shape);
}

template<
    typename Layout,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_static_extent_spec_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    return make_mdspan_impl<Layout, view_accessor_t<View>>(view, input_shape);
}

template<
    typename StaticExtentsSpec,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = view_rank_v<View>;
    using BaseExtents = extract_extents_t<InputExtents, Rank>;
    using Extents = static_extents_t<BaseExtents, StaticExtentsSpec>;
    return make_mdspan_impl<view_layout_t<View>, view_accessor_t<View>>(view, Extents{input_shape});
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec> && !make_mdspan_helpers::is_extent_like_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = view_rank_v<View>;
    using BaseExtents = extract_extents_t<InputExtents, Rank>;
    using Extents = static_extents_t<BaseExtents, StaticExtentsSpec>;
    return make_mdspan_impl<Layout, view_accessor_t<View>>(view, Extents{input_shape});
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    template <typename> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    constexpr std::size_t Rank = view_rank_v<View>;
    using BaseExtents = extract_extents_t<InputExtents, Rank>;
    using Extents = static_extents_t<BaseExtents, StaticExtentsSpec>;
    return make_mdspan_impl<Layout, Accessor>(view, Extents{input_shape});
}

template<
    typename StaticExtentsSpec,
    template <typename> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    constexpr std::size_t Rank = view_rank_v<View>;
    using BaseExtents = extract_extents_t<InputExtents, Rank>;
    using Extents = static_extents_t<BaseExtents, StaticExtentsSpec>;
    return make_mdspan_impl<view_layout_t<View>, Accessor>(view, Extents{input_shape});
}

template<
    typename StaticExtentsSpec,
    template <typename, std::size_t> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = typename bind_default_alignment<AccessorPolicy>::template type<Value>;
    constexpr std::size_t Rank = view_rank_v<View>;
    using BaseExtents = extract_extents_t<InputExtents, Rank>;
    using Extents = static_extents_t<BaseExtents, StaticExtentsSpec>;
    return make_mdspan_impl<view_layout_t<View>, Accessor>(view, Extents{input_shape});
}

template<
    typename StaticExtentsSpec,
    typename Layout,
    template <typename, std::size_t> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = typename bind_default_alignment<AccessorPolicy>::template type<Value>;
    constexpr std::size_t Rank = view_rank_v<View>;
    using BaseExtents = extract_extents_t<InputExtents, Rank>;
    using Extents = static_extents_t<BaseExtents, StaticExtentsSpec>;
    return make_mdspan_impl<Layout, Accessor>(view, Extents{input_shape});
}

template<
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    return make_mdspan_impl<view_layout_t<View>, view_accessor_t<View>>(view, input_shape);
}

template<
    template <typename> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    return make_mdspan<view_layout_t<View>, AccessorPolicy>(view, input_shape);
}

template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    return make_mdspan<view_layout_t<View>, bind_default_alignment<AccessorPolicy>::template type>(view, input_shape);
}

template<
    typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_policy_v<AccessorPolicy> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = accessor_policy_t<AccessorPolicy, Value>;
    return make_mdspan_impl<view_layout_t<View>, Accessor>(view, input_shape);
}

template<
    typename Layout,
    typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_static_extent_spec_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              make_mdspan_helpers::is_accessor_policy_v<AccessorPolicy> &&
                              !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = accessor_policy_t<AccessorPolicy, Value>;
    return make_mdspan_impl<Layout, Accessor>(view, input_shape);
}

template<
    typename Accessor,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_type_v<Accessor> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    return make_mdspan_impl<view_layout_t<View>, Accessor>(view, input_shape);
}

template<
    typename Layout,
    typename Accessor,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_static_extent_spec_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              make_mdspan_helpers::is_accessor_type_v<Accessor> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    return make_mdspan_impl<Layout, Accessor>(view, input_shape);
}

template<
    typename Layout = layout_right,
    template <typename> typename AccessorPolicy = restrict_accessor,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_container_like_v<Container> && !make_mdspan_helpers::is_extent_like_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = AccessorPolicy<Value>;
    return make_mdspan_container_impl<Extents, Layout, Accessor>(container, Extents{input_shape});
}

template<
    template <typename> typename AccessorPolicy,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    return make_mdspan<layout_right, AccessorPolicy>(container, input_shape);
}

template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    return make_mdspan<layout_right, bind_default_alignment<AccessorPolicy>::template type>(container, input_shape);
}

template<
    typename AccessorPolicy,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_policy_v<AccessorPolicy> && make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = accessor_policy_t<AccessorPolicy, Value>;
    return make_mdspan_container_impl<Extents, layout_right, Accessor>(container, Extents{input_shape});
}

template<
    typename Layout,
    typename AccessorPolicy,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_accessor_policy_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_type_v<Layout> && make_mdspan_helpers::is_accessor_policy_v<AccessorPolicy> &&
                              make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = accessor_policy_t<AccessorPolicy, Value>;
    return make_mdspan_container_impl<Extents, Layout, Accessor>(container, Extents{input_shape});
}

template<
    typename Accessor,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_type_v<Accessor> && make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    return make_mdspan_container_impl<Extents, layout_right, Accessor>(container, Extents{input_shape});
}

template<
    typename Layout,
    typename Accessor,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_accessor_policy_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_type_v<Layout> && make_mdspan_helpers::is_accessor_type_v<Accessor> &&
                              make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    return make_mdspan_container_impl<Extents, Layout, Accessor>(container, Extents{input_shape});
}

template<
    typename Layout = layout_right,
    template <typename> typename AccessorPolicy = restrict_accessor,
    typename Value,
    typename InputExtents,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_accessor_policy_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_type_v<Layout>, int> = 0
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = AccessorPolicy<Value>;
    return make_mdspan_pointer_impl<Extents, Layout, Accessor>(data, Extents{input_shape});
}

template<
    typename AccessorPolicy,
    typename Value,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_policy_v<AccessorPolicy>, int> = 0
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = accessor_policy_t<AccessorPolicy, Value>;
    return make_mdspan_pointer_impl<Extents, layout_right, Accessor>(data, Extents{input_shape});
}

template<
    typename Layout,
    typename AccessorPolicy,
    typename Value,
    typename InputExtents,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_accessor_policy_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_type_v<Layout> && make_mdspan_helpers::is_accessor_policy_v<AccessorPolicy>, int> = 0
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = accessor_policy_t<AccessorPolicy, Value>;
    return make_mdspan_pointer_impl<Extents, Layout, Accessor>(data, Extents{input_shape});
}

template<
    typename Accessor,
    typename Value,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_type_v<Accessor>, int> = 0
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    return make_mdspan_pointer_impl<Extents, layout_right, Accessor>(data, Extents{input_shape});
}

template<
    typename Layout,
    typename Accessor,
    typename Value,
    typename InputExtents,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_accessor_policy_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_type_v<Layout> && make_mdspan_helpers::is_accessor_type_v<Accessor>, int> = 0
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    return make_mdspan_pointer_impl<Extents, Layout, Accessor>(data, Extents{input_shape});
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
    using namespace make_mdspan_helpers;
    return make_mdspan<layout_right, bind_default_alignment<AccessorPolicy>::template type>(data, input_shape);
}

}  // namespace atlas
