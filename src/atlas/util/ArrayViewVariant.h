#pragma once
#pragma once

#include <string>
#include <variant>

#include "atlas/array.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace util {

namespace detail {

using namespace array;

// Container struct for a list of types.
template <typename... Types>
struct TypeList {};

// Container struct for a list of integers.
template <int... Ints>
struct IntList {};

// Check if type is a TypeList.
template <typename...>
struct IsTypeList : std::false_type {};

template <typename... Types>
struct IsTypeList<TypeList<Types...>> : std::true_type {};

// Check if type is an IntList.
template <typename...>
struct IsIntList : std::false_type {};

template <int... Ints>
struct IsIntList<IntList<Ints...>> : std::true_type {};

// Concatenate two TypeLisTypes.
template <typename... Types, typename... OtherTypes>
constexpr auto concatenateTypeLisTypes(TypeList<Types...>, TypeList<OtherTypes...>) {
  return TypeList<Types..., OtherTypes...>{};
}

// Convert types in a TypeList to const types.
template <typename... Types>
constexpr auto constTypeListImpl(TypeList<Types...>) {
  return TypeList<std::add_const_t<Types>...>{};
}

template <typename Types>
using ConstTypeList = decltype(constTypeListImpl(Types{}));

// Convert a TypeList to a std::variant of the same types.
template <typename... Types>
struct TypeListToVariantImpl {
  constexpr TypeListToVariantImpl(TypeList<Types...>) {};
  using variant = std::variant<Types...>;
};

template <typename Types>
using TypeListToVariant =
    typename decltype(TypeListToVariantImpl(Types{}))::variant;

// Iterate over each rank in an IntList. Create a TypeList of ArrayViews.
template <typename Value, int Rank, int... Ranks>
constexpr auto forEachRank(TypeList<Value>, IntList<Rank, Ranks...>) {
  return concatenateTypeLisTypes(
      TypeList<ArrayView<Value, Rank>>{},
      forEachRank(TypeList<Value>{}, IntList<Ranks...>{}));
}

template <typename Value>
constexpr auto forEachRank(TypeList<Value>, IntList<>) {
  return TypeList<>{};
}

// Iterate over each value-type in a TypeList and each rank in an IntList.
// Create a TypeList of ArrayViews.
template <typename Value, typename... Values, int... Ranks>
constexpr auto forEachValue(TypeList<Value, Values...>, IntList<Ranks...>) {
  return concatenateTypeLisTypes(
      forEachRank(TypeList<Value>{}, IntList<Ranks...>{}),
      forEachValue(TypeList<Values...>{}, IntList<Ranks...>{}));
}

template <int... Ranks>
constexpr auto forEachValue(TypeList<>, IntList<Ranks...>) {
  return TypeList<>{};
}

// Get TypeList of ArrayViews.
template <typename Values, typename Ranks>
using GetArrayViews = decltype(forEachValue(Values{}, Ranks{}));

// Match an Array to the ArrayView types in variant. Once found, return
// ArrayView. Throw exception if Array is not matched.
template <typename Variant, typename Arr, typename View, typename... Views>
Variant setVariant(Arr& arr, TypeList<View, Views...>) {
  constexpr auto Rank = View::rank();
  using Value = typename View::value_type;

  if (arr.rank() == Rank && array::DataType::kind<Value>() == arr.datatype()) {
    return Variant{array::make_view<Value, Rank>(arr)};
  }
  return setVariant<Variant>(arr, TypeList<Views...>{});
}

template <typename Variant, typename Arr>
Variant setVariant(Arr& arr, TypeList<>) {
  const auto datatypeStr = array::DataType::kind_to_str(arr.datatype().kind());
  const auto rankStr = std::to_string(arr.rank());
  throw_Exception("Could not find ArrayView<" + datatypeStr + ", " + rankStr +
                      "> in variant type.",
                  Here());
}

}  // namespace detail

/// @brief    Helper struct to contain a list of array value-types.
template <typename... Types>
using Values = detail::TypeList<Types...>;

/// @brief    Helper struct to contain a list of array ranks.
template <int... Ints>
using Ranks = detail::IntList<Ints...>;

/// @brief    Default value-types for ArrayView variant.
using DefaultValues = Values<float, double>;

/// @brief    Default ranks for ArrayView variant.
using DefaultRanks = Ranks<1, 2, 3>;


/// @ brief   Takes an array and creates a std::variant of possible ArrayViews.
///
/// @details Possible array views can be set by providing Values<types...> and
///          Ranks<inTypes...> as template argumenTypes.
template <typename ValuesList = DefaultValues,
          typename RanksList = DefaultRanks,
          typename = std::enable_if_t<detail::IsTypeList<ValuesList>::value>,
          typename = std::enable_if_t<detail::IsIntList<RanksList>::value>>
auto make_array_view_variant(array::Array& array) {
  using ArrayViewList = detail::GetArrayViews<ValuesList, RanksList>;
  using Variant = detail::TypeListToVariant<ArrayViewList>;
  return detail::setVariant<Variant>(array, ArrayViewList{});
}

/// @ brief   Takes a const array and creates a std::variant of possible
///           ArrayViews.
///
/// @details Possible array views can be set by providing Values<types...> and
///          Ranks<ints...> as template arguments. Value-types are converted to
///          const before constructing the ArrayView variant.
template <typename ValuesList = DefaultValues,
          typename RanksList = DefaultRanks,
          typename = std::enable_if_t<detail::IsTypeList<ValuesList>::value>,
          typename = std::enable_if_t<detail::IsIntList<RanksList>::value>>
auto make_array_view_variant(const array::Array& array) {
  using ConstValuesList = detail::ConstTypeList<ValuesList>;
  using ArrayViewList = detail::GetArrayViews<ConstValuesList, RanksList>;
  using Variant = detail::TypeListToVariant<ArrayViewList>;
  return detail::setVariant<Variant>(array, ArrayViewList{});
}

}  // namespace util
}  // namespace atlas
