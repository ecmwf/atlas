/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <type_traits>
#include <variant>

#include "atlas/array/ArrayView.h"

namespace atlas {
namespace array {

namespace detail {

// Container struct for a list of types.
template <typename... Ts>
struct Types {
  using add_const = Types<std::add_const_t<Ts>...>;
};

// Container struct for a list of integers.
template <int... Is>
struct Ints {};

template <typename ValueTypes, typename Ranks, typename... ArrayViews>
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

template <typename ValueTypes, typename Ranks>
using Variant = typename VariantHelper<ValueTypes, Ranks>::type;

using VariantValueTypes =
    detail::Types<float, double, int, long, unsigned long>;

using VariantRanks = detail::Ints<1, 2, 3, 4, 5, 6, 7, 8, 9>;

}  // namespace detail

class Array;

/// @brief Variant containing all supported non-const ArrayView alternatives.
using ArrayViewVariant =
    detail::Variant<detail::VariantValueTypes, detail::VariantRanks>;

/// @brief Variant containing all supported const ArrayView alternatives.
using ConstArrayViewVariant =
    detail::Variant<detail::VariantValueTypes::add_const, detail::VariantRanks>;

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
template <int... Ranks, typename View>
constexpr bool is_rank(const View&) {
  return ((std::decay_t<View>::rank() == Ranks) || ...);
}
/// @brief Return true if View::value_type is any of ValuesTypes...
template <typename... ValueTypes, typename View>
constexpr bool is_value_type(const View&) {
  using ValueType = typename std::decay_t<View>::value_type;
  return ((std::is_same_v<ValueType, ValueTypes>) || ...);
}

/// @brief Return true if View::non_const_value_type is any of ValuesTypes...
template <typename... ValueTypes, typename View>
constexpr bool is_non_const_value_type(const View&) {
  using ValueType = typename std::decay_t<View>::non_const_value_type;
  return ((std::is_same_v<ValueType, ValueTypes>) || ...);
}

}  // namespace array
}  // namespace atlas
