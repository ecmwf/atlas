/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <variant>

#include "atlas/array.h"

namespace atlas {
namespace array {

namespace detail {

using namespace array;

// Container struct for a list of types.
template <typename... Ts>
struct Types {};

// Container struct for a list of integers.
template <int... Is>
struct Ints {};

template <typename...>
struct VariantHelper;

// Recursively construct ArrayView std::variant from types Ts and ranks Is.
template <typename... ArrayViews, typename T, typename... Ts, int... Is>
struct VariantHelper<Types<ArrayViews...>, Types<T, Ts...>, Ints<Is...>> {
  using type = typename VariantHelper<
      Types<ArrayViews..., ArrayView<const T, Is>..., ArrayView<T, Is>...>,
      Types<Ts...>, Ints<Is...>>::type;
};

// End recursion.
template <typename... ArrayViews, int... Is>
struct VariantHelper<Types<ArrayViews...>, Types<>, Ints<Is...>> {
  using type = std::variant<ArrayViews...>;
};

template <typename Values, typename Ranks>
using Variant = typename VariantHelper<Types<>, Values, Ranks>::type;

}  // namespace detail

/// @brief Supported ArrayView value types.
using Values = detail::Types<float, double, int, long, unsigned long>;

/// @brief Supported ArrayView ranks.
using Ranks = detail::Ints<1, 2, 3, 4, 5, 6, 7, 8, 9>;

/// @brief Variant containing all supported ArrayView alternatives.
using ArrayViewVariant = detail::Variant<Values, Ranks>;

/// @brief Create an ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_view_variant(Array& array);

/// @brief Create a const ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_view_variant(const Array& array);

/// @brief Create a host ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_host_view_variant(Array& array);

/// @brief Create a host const ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_host_view_variant(const Array& array);

/// @brief Create a device ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_device_view_variant(Array& array);

/// @brief Create a const device ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_device_view_variant(const Array& array);

}  // namespace array
}  // namespace atlas
