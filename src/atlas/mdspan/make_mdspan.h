#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

#include "atlas/mdspan/make_mdspan_helpers.h"
#include "atlas/mdspan/mdspan.h"
#include "atlas/mdspan/mdspan_introspection.h"

namespace atlas {

//------------------------------------------------------------------------------
/// @name make_mdspan
///
/// @brief Adapt Atlas views, mdspans, containers, and raw pointers to mdspan.
///
/// @details `make_mdspan` is the public conversion and reshaping entry point for
/// mdspan-compatible storage in Atlas. It creates a non-owning mdspan over the
/// input data and derives as much metadata as possible from the input object: view
/// overloads preserve the view data handle, extents, strides, layout, and accessor;
/// container and raw-pointer overloads assume contiguous storage and use
/// layout_right unless another layout is requested.
///
/// The overload set is organised around four common decisions:
///
/// - **Input source**: pass an Atlas view, an existing mdspan, a contiguous
///   container with `data()` and `size()`, or a raw pointer.
/// - **Shape**: omit the shape when it can be read from the input, pass a runtime
///   shape such as `std::array<std::size_t, Rank>` or `dims<Rank>`, or provide an
///   explicit mdspan extents type such as `dims<2>` or
///   `extents<std::size_t, dynamic_extent, 3>`.
/// - **Layout**: let the input layout be preserved for views/mdspans, use the
///   contiguous layout_right default for containers and pointers, or specify a
///   layout such as `layout_stride` or `layout_right` explicitly.
/// - **Accessor**: use the default restrict accessor, choose an accessor template
///   such as `default_accessor` or `aligned_accessor`, pass an accessor policy
///   such as `aligned_accessor_policy<64>`, or pass a fully bound accessor type
///   such as `aligned_accessor<double, 64>`.
///
/// Static extent helpers refine the resulting mdspan type while keeping the call
/// site compact. For example `extent_with_static_last_dim<3>` keeps the last dimension fixed
/// at compile time, while dynamic dimensions are still taken from the input shape.
///
/// @code
/// auto view_span      = make_mdspan(view);
/// auto layout_span    = make_mdspan<layout_right>(view);
/// auto static_span    = make_mdspan<extent_with_static_last_dim<3>>(view);
/// auto accessor_span  = make_mdspan<aligned_accessor_policy<64>>(view);
/// auto vector_span    = make_mdspan(values, std::array<std::size_t, 2>{2, 3});
/// auto pointer_span   = make_mdspan(data, dims<2>{2, 3});
/// auto explicit_shape = make_mdspan<dims<2>>(data, {2, 3});
/// @endcode
///
/// A typical SIMD-oriented use is to expose a fast path for a blocked field whose
/// last dimension is 32. The checks make the assumptions explicit before
/// choosing layout_right, a fixed last extent of 32, and a 64-byte
/// restrict_aligned accessor. With those guarantees, loops over the last
/// dimension have a fixed trip count and aligned non-aliased accesses.
///
/// @code
/// if (can_use_layout<layout_right>(view) && is_last_dimension_aligned(view, 64)) {
///     using layout = layout_right;
///     using accessor_policy = restrict_aligned_accessor_policy<64>;
///     switch(last_extent(view)) {
///         case 16:
///             return kernel(make_mdspan<extent_with_static_last_dim<16>, layout, accessor_policy>(view));
///         case 32:
///             return kernel(make_mdspan<extent_with_static_last_dim<32>, layout, accessor_policy>(view));
///         case 64:
///             return kernel(make_mdspan<extent_with_static_last_dim<64>, layout, accessor_policy>(view));
///         default:
///             return kernel(make_mdspan<layout, accessor_policy>(view));
///     }
/// }
/// return kernel(make_mdspan(view));
/// @endcode
///
/// The returned mdspan does not own the underlying memory. The caller must ensure
/// that the input storage outlives the mdspan and that any declared alignment,
/// restrict, layout, or static-extent assumptions are valid for the data.
///
/// @{

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with explicit extents and an optional layout/accessor policy.
///
/// @details Motivation: use this overload when the extents type is already known
/// and should be part of the resulting mdspan type, while preserving the view data
/// handle and deriving the accessor element type from the view.
///
/// @code
/// auto span = make_mdspan<extents<std::size_t, 2, 3>, layout_stride>(view);
/// @endcode
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

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with explicit extents and a two-argument accessor template.
///
/// @details Motivation: support accessors such as aligned_accessor<T, N> without
/// requiring callers to spell the element type or default alignment. The accessor
/// alignment defaults to alignof(value_type).
///
/// @code
/// auto span = make_mdspan<extents<std::size_t, 2, 3>, layout_stride, aligned_accessor>(view);
/// @endcode
template<
    typename Extents,
    typename Layout,
    template <typename, std::size_t> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_extent_like_v<Extents> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = AccessorPolicy<Value, default_accessor_alignment_v<Value>>;
    return make_mdspan_impl<Extents, Layout, Accessor>(view);
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with explicit extents and a bound accessor or accessor policy.
///
/// @details Motivation: use this overload when the accessor has already been
/// selected, for example a fully bound aligned accessor or an accessor policy
/// that carries compile-time configuration.
///
/// @code
/// auto span = make_mdspan<extents<std::size_t, 2, 3>, layout_right,
///                         aligned_accessor_policy<64>>(view);
/// @endcode
template<
    typename Extents,
    typename Layout,
    typename AccessorOrPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_extent_like_v<Extents> && make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy> &&
                              !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = accessor_or_policy_t<AccessorOrPolicy, Value>;
    return make_mdspan_impl<Extents, Layout, Accessor>(view);
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with explicit extents and an optional layout/accessor policy.
///
/// @details Motivation: let callers override the runtime shape while fixing the
/// mdspan extents type. This is useful when the view storage is valid but the
/// logical shape should be supplied by the caller.
///
/// @code
/// auto span = make_mdspan<dextents<std::size_t, 2>, layout_right>(view, shape);
/// @endcode
template<
    typename Extents,
    typename Layout = layout_stride,
    template <typename> typename AccessorPolicy = restrict_accessor,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_extent_like_v<Extents> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = AccessorPolicy<Value>;
    return make_mdspan_impl<Layout, Accessor>(view, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with explicit extents and a two-argument accessor template.
///
/// @details Motivation: combine caller-provided shape with default-aligned
/// accessor templates such as aligned_accessor<T, N>, without introducing an
/// adapter type in the public template arguments.
///
/// @code
/// auto span = make_mdspan<dextents<std::size_t, 2>, layout_stride,
///                         aligned_accessor>(view, shape);
/// @endcode
template<
    typename Extents,
    typename Layout,
    template <typename, std::size_t> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_extent_like_v<Extents> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = AccessorPolicy<Value, default_accessor_alignment_v<Value>>;
    return make_mdspan_impl<Layout, Accessor>(view, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with explicit extents and a bound accessor or accessor policy.
///
/// @details Motivation: keep the shaped extents path available when the accessor
/// is already a concrete type or a policy wrapper with compile-time parameters.
///
/// @code
/// auto span = make_mdspan<dextents<std::size_t, 2>, layout_stride,
///                         restrict_aligned_accessor<double, 64>>(view, shape);
/// @endcode
template<
    typename Extents,
    typename Layout,
    typename AccessorOrPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_extent_like_v<Extents> && make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy> &&
                              !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = accessor_or_policy_t<AccessorOrPolicy, Value>;
    return make_mdspan_impl<Layout, Accessor>(view, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view using one mdspan option.
///
/// @details Motivation: provide a compact syntax for the common cases where the
/// caller wants to specify only one of layout, static extent spec, accessor type,
/// or accessor policy, and let the view provide the remaining choices.
///
/// @code
/// auto layout_span = make_mdspan<layout_right>(view);
/// auto static_span = make_mdspan<extent_with_static_last_dim<3>>(view);
/// auto accessor_span = make_mdspan<aligned_accessor<double, 64>>(view);
/// @endcode
template<
    typename Option,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Option> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    if constexpr (is_static_extent_spec_v<Option>) {
        using Extents = static_extents_t<view_extents_t<View>, Option>;
        return make_mdspan<Extents, view_layout_t<View>, view_accessor_t<View>>(view);
    }
    else if constexpr (is_accessor_or_policy_v<Option>) {
        using Value = view_value_t<View>;
        using Accessor = accessor_or_policy_t<Option, Value>;
        return make_mdspan_impl<view_extents_t<View>, view_layout_t<View>, Accessor>(view);
    }
    else {
        return make_mdspan_impl<view_extents_t<View>, Option, view_accessor_t<View>>(view);
    }
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with explicit layout and one-argument accessor policy.
///
/// @details Motivation: select both layout and accessor family while still
/// deriving extents and element type from the view.
///
/// @code
/// auto span = make_mdspan<layout_right, default_accessor>(view);
/// @endcode
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

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with a static extent spec and explicit layout.
///
/// @details Motivation: make selected dimensions static in the mdspan type while
/// preserving the view accessor and choosing a layout explicitly.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<3>, layout_stride>(view);
/// @endcode
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
    return make_mdspan<Extents, Layout, view_accessor_t<View>>(view);
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with static extents, explicit layout, and one-argument accessor policy.
///
/// @details Motivation: combine the compile-time shape refinement of
/// StaticExtent/extent_with_static_last_dim with a chosen layout and accessor family.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<3>, layout_right, default_accessor>(view);
/// @endcode
template<
    typename StaticExtentsSpec,
    typename Layout,
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan<Extents, Layout, AccessorPolicy>(view);
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with static extents, explicit layout, and bound accessor or policy.
///
/// @details Motivation: combine static extent refinement with a layout override
/// and an accessor policy carrying compile-time configuration, such as alignment.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<32>, layout_right,
///                         restrict_aligned_accessor_policy<64>>(view);
/// @endcode
template<
    typename StaticExtentsSpec,
    typename Layout,
    typename AccessorOrPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec> &&
                              make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan<Extents, Layout, AccessorOrPolicy>(view);
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with static extents and one-argument accessor policy.
///
/// @details Motivation: refine the extents type while keeping the view layout and
/// selecting an accessor family such as default_accessor or restrict_accessor.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<3>, default_accessor>(view);
/// @endcode
template<
    typename StaticExtentsSpec,
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan<Extents, view_layout_t<View>, AccessorPolicy>(view);
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with static extents and a two-argument accessor template.
///
/// @details Motivation: use aligned-style accessor templates together with static
/// extent refinement. The default alignment is selected from the view value type.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<3>, aligned_accessor>(view);
/// @endcode
template<
    typename StaticExtentsSpec,
    template <typename, std::size_t> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan<Extents, view_layout_t<View>, AccessorPolicy>(view);
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with static extents, explicit layout, and a two-argument accessor template.
///
/// @details Motivation: give full control over static extents, layout, and
/// aligned-style accessor family while still deriving element type from the view.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<3>, layout_stride, aligned_accessor>(view);
/// @endcode
template<
    typename StaticExtentsSpec,
    typename Layout,
    template <typename, std::size_t> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Extents = static_extents_t<view_extents_t<View>, StaticExtentsSpec>;
    return make_mdspan<Extents, Layout, AccessorPolicy>(view);
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with a one-argument accessor policy.
///
/// @details Motivation: provide accessor-first syntax for the common case where
/// extents and layout should come from the view.
///
/// @code
/// auto span = make_mdspan<default_accessor>(view);
/// @endcode
template<
    template <typename> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    return make_mdspan<view_extents_t<View>, view_layout_t<View>, AccessorPolicy>(view);
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with a two-argument accessor template.
///
/// @details Motivation: provide accessor-first syntax for aligned-style accessors
/// while deriving the element type, extents, layout, and default alignment.
///
/// @code
/// auto span = make_mdspan<aligned_accessor>(view);
/// @endcode
template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    return make_mdspan<view_extents_t<View>, view_layout_t<View>, AccessorPolicy>(view);
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view with explicit layout and bound accessor or accessor policy.
///
/// @details Motivation: use a fully specified accessor type or policy wrapper
/// together with a layout override, while deriving extents from the view.
///
/// @code
/// auto span = make_mdspan<layout_stride, aligned_accessor<double, 64>>(view);
/// @endcode
template<
    typename Layout,
    typename AccessorOrPolicy,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_static_extent_spec_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy> &&
                              !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = accessor_or_policy_t<AccessorOrPolicy, Value>;
    return make_mdspan_impl<view_extents_t<View>, Layout, Accessor>(view);
}

//------------------------------------------------------------------------------
/// @brief Create an mdspan from a view using the view extents, layout, and accessor.
///
/// @details Motivation: this is the shortest conversion from an Atlas view-like
/// object to mdspan when the view already carries the desired metadata.
///
/// @code
/// auto span = make_mdspan(view);
/// @endcode
template<typename View, typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0>
auto make_mdspan(View& view) {
    using namespace make_mdspan_helpers;
    return make_mdspan_impl<view_extents_t<View>, view_layout_t<View>,
                                    view_accessor_t<View>>(view);
}

//------------------------------------------------------------------------------
/// @brief Create a rank-1 mdspan over a container with optional layout/accessor policy.
///
/// @details Motivation: adapt contiguous container-like objects, such as
/// std::vector, to mdspan without requiring the caller to provide a shape.
///
/// @code
/// std::vector<double> values(6);
/// auto span = make_mdspan(values);
/// @endcode
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

//------------------------------------------------------------------------------
/// @brief Create a rank-2 mdspan over a vector of fixed-size arrays.
///
/// @details Motivation: adapt common connectivity storage such as
/// `std::vector<std::array<T, N>>` while preserving the fixed inner extent in the
/// mdspan type.
///
/// @code
/// std::vector<std::array<gidx_t, 3>> cells;
/// auto span = make_mdspan(cells);
/// @endcode
template <typename T, std::size_t N>
auto make_mdspan(std::vector<std::array<T, N>>& vector) {
    return mdspan<T, extents<std::size_t, dynamic_extent, N>>{reinterpret_cast<T*>(vector.data()), vector.size()};
}

//------------------------------------------------------------------------------
/// @brief Create a rank-1 mdspan over a container with a one-argument accessor policy.
///
/// @details Motivation: keep the convenient container syntax while letting the
/// caller choose an accessor family.
///
/// @code
/// auto span = make_mdspan<default_accessor>(values);
/// @endcode
template<
    template <typename> typename AccessorPolicy,
    typename Container,
    typename std::enable_if_t<make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    return make_mdspan<layout_right, AccessorPolicy>(container);
}

//------------------------------------------------------------------------------
/// @brief Create a rank-1 mdspan over a container with a two-argument accessor template.
///
/// @details Motivation: use aligned-style accessor templates for containers while
/// deriving the element type and default alignment from container.data().
///
/// @code
/// auto span = make_mdspan<aligned_accessor>(values);
/// @endcode
template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename Container,
    typename std::enable_if_t<make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    using Extents = dextents<std::size_t, 1>;
    using Accessor = AccessorPolicy<Value, default_accessor_alignment_v<Value>>;
    return make_mdspan_container_impl<Extents, layout_right, Accessor>(container, Extents{container.size()});
}

//------------------------------------------------------------------------------
/// @brief Create a rank-1 mdspan over a container with a bound accessor or accessor policy.
///
/// @details Motivation: allow concrete accessor types and policy wrappers with
/// explicit compile-time parameters for container-backed mdspans.
///
/// @code
/// auto span = make_mdspan<aligned_accessor<double, 64>>(values);
/// @endcode
template<
    typename AccessorOrPolicy,
    typename Container,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy> &&
                              make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    using Extents = dextents<std::size_t, 1>;
    using Accessor = accessor_or_policy_t<AccessorOrPolicy, Value>;
    return make_mdspan_container_impl<Extents, layout_right, Accessor>(container, Extents{container.size()});
}

//------------------------------------------------------------------------------
/// @brief Create a rank-1 mdspan over a container with explicit layout and bound accessor or policy.
///
/// @details Motivation: combine a layout override with a concrete accessor choice
/// for contiguous container data.
///
/// @code
/// auto span = make_mdspan<layout_right, default_accessor<double>>(values);
/// @endcode
template<
    typename Layout,
    typename AccessorOrPolicy,
    typename Container,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_accessor_policy_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy> &&
                              make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    using Extents = dextents<std::size_t, 1>;
    using Accessor = accessor_or_policy_t<AccessorOrPolicy, Value>;
    return make_mdspan_container_impl<Extents, Layout, Accessor>(container, Extents{container.size()});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with explicit layout and one-argument accessor policy.
///
/// @details Motivation: let the caller supply a logical shape while selecting
/// layout and accessor family explicitly.
///
/// @code
/// auto span = make_mdspan<layout_right, default_accessor>(view, shape);
/// @endcode
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

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view using one mdspan option.
///
/// @details Motivation: provide compact shaped syntax when the caller only needs
/// to specify one option: layout, static extent spec, accessor type, or policy.
///
/// @code
/// auto layout_span = make_mdspan<layout_right>(view, shape);
/// auto static_span = make_mdspan<extent_with_static_last_dim<3>>(view, shape);
/// @endcode
template<
    typename Option,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Option> && !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    if constexpr (is_static_extent_spec_v<Option>) {
        constexpr std::size_t Rank = view_rank_v<View>;
        using BaseExtents = extract_extents_t<InputExtents, Rank>;
        using Extents = static_extents_t<BaseExtents, Option>;
        return make_mdspan<Extents, view_layout_t<View>, view_accessor_t<View>>(view, Extents{input_shape});
    }
    else if constexpr (is_accessor_or_policy_v<Option>) {
        using Value = view_value_t<View>;
        using Accessor = accessor_or_policy_t<Option, Value>;
        return make_mdspan_impl<view_layout_t<View>, Accessor>(view, input_shape);
    }
    else {
        return make_mdspan_impl<Option, view_accessor_t<View>>(view, input_shape);
    }
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with static extents and explicit layout.
///
/// @details Motivation: convert selected runtime extents to static mdspan extents
/// while using a caller-chosen layout.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<3>, layout_stride>(view, shape);
/// @endcode
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
    return make_mdspan<Extents, Layout, view_accessor_t<View>>(view, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with static extents, layout, and one-argument accessor policy.
///
/// @details Motivation: combine caller-provided shape, static extent refinement,
/// layout selection, and accessor-family selection.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<3>, layout_right, default_accessor>(view, shape);
/// @endcode
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
    constexpr std::size_t Rank = view_rank_v<View>;
    using BaseExtents = extract_extents_t<InputExtents, Rank>;
    using Extents = static_extents_t<BaseExtents, StaticExtentsSpec>;
    return make_mdspan<Extents, Layout, AccessorPolicy>(view, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with static extents and one-argument accessor policy.
///
/// @details Motivation: keep the view layout while refining the extents type and
/// selecting an accessor family.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<3>, default_accessor>(view, shape);
/// @endcode
template<
    typename StaticExtentsSpec,
    template <typename> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = view_rank_v<View>;
    using BaseExtents = extract_extents_t<InputExtents, Rank>;
    using Extents = static_extents_t<BaseExtents, StaticExtentsSpec>;
    return make_mdspan<Extents, view_layout_t<View>, AccessorPolicy>(view, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with static extents and a two-argument accessor template.
///
/// @details Motivation: use aligned-style accessor templates on shaped views
/// while deriving element type and default alignment from the view.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<3>, aligned_accessor>(view, shape);
/// @endcode
template<
    typename StaticExtentsSpec,
    template <typename, std::size_t> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<make_mdspan_helpers::is_static_extent_spec_v<StaticExtentsSpec>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = view_rank_v<View>;
    using BaseExtents = extract_extents_t<InputExtents, Rank>;
    using Extents = static_extents_t<BaseExtents, StaticExtentsSpec>;
    return make_mdspan<Extents, view_layout_t<View>, AccessorPolicy>(view, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with static extents, layout, and two-argument accessor template.
///
/// @details Motivation: provide the full shaped static-extents syntax for
/// aligned-style accessor templates.
///
/// @code
/// auto span = make_mdspan<extent_with_static_last_dim<3>, layout_stride, aligned_accessor>(view, shape);
/// @endcode
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
    constexpr std::size_t Rank = view_rank_v<View>;
    using BaseExtents = extract_extents_t<InputExtents, Rank>;
    using Extents = static_extents_t<BaseExtents, StaticExtentsSpec>;
    return make_mdspan<Extents, Layout, AccessorPolicy>(view, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view using the view layout and accessor.
///
/// @details Motivation: reinterpret a view with a caller-supplied logical shape
/// without changing layout or accessor choices.
///
/// @code
/// auto span = make_mdspan(view, shape);
/// @endcode
template<
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    return make_mdspan_impl<view_layout_t<View>, view_accessor_t<View>>(view, input_shape);
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with a one-argument accessor policy.
///
/// @details Motivation: provide accessor-first shaped syntax while keeping the
/// view layout and deriving the extents from the supplied shape.
///
/// @code
/// auto span = make_mdspan<default_accessor>(view, shape);
/// @endcode
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

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with a two-argument accessor template.
///
/// @details Motivation: provide accessor-first shaped syntax for aligned-style
/// accessors using the view value type and default alignment.
///
/// @code
/// auto span = make_mdspan<aligned_accessor>(view, shape);
/// @endcode
template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    return make_mdspan<view_extents_t<View>, view_layout_t<View>, AccessorPolicy>(view, input_shape);
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan from a view with explicit layout and bound accessor or policy.
///
/// @details Motivation: allow the shaped view path to use a concrete accessor or
/// policy wrapper together with a layout override.
///
/// @code
/// auto span = make_mdspan<layout_stride, aligned_accessor<double, 64>>(view, shape);
/// @endcode
template<
    typename Layout,
    typename AccessorOrPolicy,
    typename InputExtents,
    typename View,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_static_extent_spec_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy> &&
                              !make_mdspan_helpers::is_container_like_v<View>, int> = 0
>
auto make_mdspan(View& view, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = view_value_t<View>;
    using Accessor = accessor_or_policy_t<AccessorOrPolicy, Value>;
    return make_mdspan_impl<Layout, Accessor>(view, input_shape);
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a container with optional layout/accessor policy.
///
/// @details Motivation: adapt contiguous container storage to an arbitrary logical
/// shape supplied by the caller.
///
/// @code
/// auto span = make_mdspan(values, std::array<std::size_t, 2>{2, 3});
/// @endcode
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

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a container with a one-argument accessor policy.
///
/// @details Motivation: keep shaped container syntax concise while selecting an
/// accessor family.
///
/// @code
/// auto span = make_mdspan<default_accessor>(values, shape);
/// @endcode
template<
    template <typename> typename AccessorPolicy,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    return make_mdspan<layout_right, AccessorPolicy>(container, input_shape);
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a container with a two-argument accessor template.
///
/// @details Motivation: use aligned-style accessors for shaped container data,
/// with element type and default alignment derived from the container.
///
/// @code
/// auto span = make_mdspan<aligned_accessor>(values, shape);
/// @endcode
template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = AccessorPolicy<Value, default_accessor_alignment_v<Value>>;
    return make_mdspan_container_impl<Extents, layout_right, Accessor>(container, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a container with a bound accessor or accessor policy.
///
/// @details Motivation: allow explicit accessor types and policy wrappers for
/// shaped container-backed mdspans.
///
/// @code
/// auto span = make_mdspan<aligned_accessor<double, 64>>(values, shape);
/// @endcode
template<
    typename AccessorOrPolicy,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy> &&
                              make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = accessor_or_policy_t<AccessorOrPolicy, Value>;
    return make_mdspan_container_impl<Extents, layout_right, Accessor>(container, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a container with explicit layout and bound accessor or policy.
///
/// @details Motivation: provide full layout and accessor control for shaped
/// container data.
///
/// @code
/// auto span = make_mdspan<layout_right, default_accessor<double>>(values, shape);
/// @endcode
template<
    typename Layout,
    typename AccessorOrPolicy,
    typename Container,
    typename InputExtents,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_accessor_policy_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy> &&
                              make_mdspan_helpers::is_container_like_v<Container>, int> = 0
>
auto make_mdspan(Container& container, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    using Value = std::remove_pointer_t<std::remove_reference_t<decltype(container.data())>>;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = accessor_or_policy_t<AccessorOrPolicy, Value>;
    return make_mdspan_container_impl<Extents, Layout, Accessor>(container, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a raw pointer with explicit extents and optional layout/accessor policy.
///
/// @details Motivation: give the runtime shape argument a concrete rank and index
/// type from the mdspan extents template argument. This supports braced shape
/// syntax while keeping the resulting mdspan extents type explicit.
///
/// @code
/// double data[6]{};
/// auto span = make_mdspan<dims<2>>(data, {2, 3});
/// @endcode
template<
    typename Extents,
    typename Layout = layout_right,
    template <typename> typename AccessorPolicy = restrict_accessor,
    typename Value,
    typename std::enable_if_t<make_mdspan_helpers::is_extent_like_v<Extents> && !make_mdspan_helpers::is_extent_like_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_policy_v<Layout> && !make_mdspan_helpers::is_accessor_type_v<Layout>, int> = 0
>
auto make_mdspan(Value* data, std::array<typename Extents::index_type, Extents::rank()> input_shape) {
    using Accessor = AccessorPolicy<Value>;
    return make_mdspan_helpers::make_mdspan_pointer_impl<Extents, Layout, Accessor>(data, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a raw pointer with optional layout/accessor policy.
///
/// @details Motivation: adapt raw contiguous storage to mdspan when no container
/// or view wrapper is available.
///
/// @code
/// double data[6]{};
/// auto span = make_mdspan(data, dims<2>{2, 3});
/// @endcode
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

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a raw pointer with a bound accessor or accessor policy.
///
/// @details Motivation: use explicit accessor types or policy wrappers with raw
/// storage, for example to express alignment or restrict semantics.
///
/// @code
/// auto span = make_mdspan<aligned_accessor<double, 64>>(data, shape);
/// @endcode
template<
    typename AccessorOrPolicy,
    typename Value,
    typename InputExtents,
    typename std::enable_if_t<make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy>, int> = 0
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = accessor_or_policy_t<AccessorOrPolicy, Value>;
    return make_mdspan_pointer_impl<Extents, layout_right, Accessor>(data, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a raw pointer with explicit layout and bound accessor or policy.
///
/// @details Motivation: provide full layout and accessor control when adapting raw
/// pointer storage.
///
/// @code
/// auto span = make_mdspan<layout_stride, default_accessor<double>>(data, shape);
/// @endcode
template<
    typename Layout,
    typename AccessorOrPolicy,
    typename Value,
    typename InputExtents,
    typename std::enable_if_t<!make_mdspan_helpers::is_extent_like_v<Layout> && !make_mdspan_helpers::is_accessor_policy_v<Layout> &&
                              !make_mdspan_helpers::is_accessor_type_v<Layout> &&
                              make_mdspan_helpers::is_accessor_or_policy_v<AccessorOrPolicy>, int> = 0
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = accessor_or_policy_t<AccessorOrPolicy, Value>;
    return make_mdspan_pointer_impl<Extents, Layout, Accessor>(data, Extents{input_shape});
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a raw pointer with a one-argument accessor policy.
///
/// @details Motivation: keep raw-pointer syntax concise while selecting an
/// accessor family and defaulting to layout_right.
///
/// @code
/// auto span = make_mdspan<default_accessor>(data, shape);
/// @endcode
template<
    template <typename> typename AccessorPolicy,
    typename Value,
    typename InputExtents
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    return make_mdspan<layout_right, AccessorPolicy>(data, input_shape);
}

//------------------------------------------------------------------------------
/// @brief Create a shaped mdspan over a raw pointer with a two-argument accessor template.
///
/// @details Motivation: use aligned-style accessor templates with raw storage,
/// deriving the element type from the pointer and the default alignment from it.
///
/// @code
/// auto span = make_mdspan<aligned_accessor>(data, shape);
/// @endcode
template<
    template <typename, std::size_t> typename AccessorPolicy,
    typename Value,
    typename InputExtents
>
auto make_mdspan(Value* data, InputExtents input_shape) {
    using namespace make_mdspan_helpers;
    constexpr std::size_t Rank = extents_rank_v<InputExtents>;
    using Extents = extract_extents_t<InputExtents, Rank>;
    using Accessor = AccessorPolicy<Value, default_accessor_alignment_v<Value>>;
    return make_mdspan_pointer_impl<Extents, layout_right, Accessor>(data, Extents{input_shape});
}

/// @}

}  // namespace atlas
