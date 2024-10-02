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
struct Types {
  using add_const = Types<std::add_const_t<Ts>...>;
};

// Container struct for a list of integers.
template <int... Is>
struct Ints {};

template <typename Values, typename Ranks, typename... ArrayViews>
struct VariantHelper;

// Recursively construct ArrayView std::variant from types Ts and ranks Is.
template <typename T, typename... Ts, int... Is, typename... ArrayViews>
struct VariantHelper<Types<T, Ts...>, Ints<Is...>, ArrayViews...> {
  using type = typename VariantHelper<Types<Ts...>, Ints<Is...>, ArrayViews...,
                                      ArrayView<T, Is>...>::type;
};

// End recursion.
template <int... Is, typename... ArrayViews>
struct VariantHelper<Types<>, Ints<Is...>, ArrayViews...> {
  using type = std::variant<ArrayViews...>;
};

template <typename Values, typename Ranks>
using Variant = typename VariantHelper<Values, Ranks>::type;

}  // namespace detail

/// @brief Supported ArrayView value types.
using ValueTypes = detail::Types<float, double, int, long, unsigned long>;

/// @brief Supported ArrayView ranks.
using Ranks = detail::Ints<1, 2, 3, 4, 5, 6, 7, 8, 9>;

/// @brief Variant containing all supported non-const ArrayView alternatives.
using ArrayViewVariant = detail::Variant<ValueTypes, Ranks>;

/// @brief Variant containing all supported const ArrayView alternatives.
using ConstArrayViewVariant = detail::Variant<ValueTypes::add_const, Ranks>;

/// @brief Create an ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_view_variant(Array& array);

/// @brief Create a const ArrayView and assign to an ArrayViewVariant.
ConstArrayViewVariant make_view_variant(const Array& array);

/// @brief Create a host ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_host_view_variant(Array& array);

/// @brief Create a const host ArrayView and assign to an ArrayViewVariant.
ConstArrayViewVariant make_host_view_variant(const Array& array);

/// @brief Create a device ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_device_view_variant(Array& array);

/// @brief Create a const device ArrayView and assign to an ArrayViewVariant.
ConstArrayViewVariant make_device_view_variant(const Array& array);

/// @brief Return true if View::rank() is any of Ranks...
template <typename View, int... Ranks>
constexpr bool is_rank() {
  return ((std::decay_t<View>::rank() == Ranks) || ...);
}

/// @brief Return true if View::value_type is any of ValuesTypes...
template <typename View, typename... ValueTypes>
constexpr bool is_value_type() {
  using ValueType = typename std::decay_t<View>::value_type;
  return ((std::is_same_v<ValueType, ValueTypes>) || ...);
}

/// @brief Return true if View::non_const_value_type is any of ValuesTypes...
template <typename View, typename... ValueTypes>
constexpr bool is_non_const_value_type() {
  using ValueType = typename std::decay_t<View>::non_const_value_type;
  return ((std::is_same_v<ValueType, ValueTypes>) || ...);
}

}  // namespace array
}  // namespace atlas
