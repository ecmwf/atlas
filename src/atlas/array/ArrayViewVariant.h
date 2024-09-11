/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <variant>

#include "atlas/array.h"
#include "eckit/utils/Overloaded.h"

namespace atlas {
namespace array {

namespace detail {

using namespace array;

// Container struct for a list of types.
template <typename... Types>
struct TypeList {};

// Container struct for a list of integers.
template <idx_t... Ints>
struct IntList {};

// Supported ArrayView value types.
constexpr auto ValueList = TypeList<float, double, int, long, unsigned long>{};

// Supported ArrayView ranks.
constexpr auto RankList = IntList<1, 2, 3, 4, 5, 6, 7, 8, 9>{};

// Helper struct to build an ArrayView variant from a list of value types and
// and a list of ranks.
template <typename... Views>
struct VariantBuilder {
  using type = std::variant<Views...>;

  // Make a VariantBuilder struct with a fully populated Views... argument.
  template <typename Value, typename... Values, idx_t... Ranks>
  static constexpr auto make(TypeList<Value, Values...>, IntList<Ranks...>) {
    return VariantBuilder<
        Views..., ArrayView<const Value, Ranks>...,
        ArrayView<Value, Ranks>...>::make(TypeList<Values...>{},
                                          IntList<Ranks...>{});
  }

  // End recursion.
  template <idx_t... Ranks>
  static constexpr VariantBuilder<Views...> make(TypeList<>,
                                                 IntList<Ranks...>) {
    return VariantBuilder<Views...>{};
  }
};

}  // namespace detail

/// @brief Variant containing all supported ArrayView alternatives.
using ArrayViewVariant = decltype(detail::VariantBuilder<>::make(
    detail::ValueList, detail::RankList))::type;

/// @brief Use overloaded pattern as visitor.
using eckit::Overloaded;

/// @brief Create an ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_view_variant(array::Array& array);

/// @brief Create a const ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_view_variant(const array::Array& array);

/// @brief Create a host ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_host_view_variant(array::Array& array);

/// @brief Create a host const ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_host_view_variant(const array::Array& array);

/// @brief Create a device ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_devive_view_variant(array::Array& array);

/// @brief Create a const device ArrayView and assign to an ArrayViewVariant.
ArrayViewVariant make_device_view_variant(const array::Array& array);

}  // namespace array
}  // namespace atlas
