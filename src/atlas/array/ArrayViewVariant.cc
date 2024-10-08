/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/ArrayViewVariant.h"

#include <string>
#include <type_traits>

#include "atlas/runtime/Exception.h"

namespace atlas {
namespace array {

using namespace detail;

namespace {

template <bool IsConst>
struct VariantTypeHelper {
  using type = ArrayViewVariant;
};

template <>
struct VariantTypeHelper<true> {
  using type = ConstArrayViewVariant;
};

template <typename ArrayType>
using VariantType =
    typename VariantTypeHelper<std::is_const_v<ArrayType>>::type;

// Match array.rank() and array.datatype() to variant types. Return result of
// makeView on a successful pattern match.
template <size_t TypeIndex = 0, typename ArrayType, typename MakeView>
VariantType<ArrayType> executeMakeView(ArrayType& array,
                                       const MakeView& makeView) {
  using View = std::variant_alternative_t<TypeIndex, VariantType<ArrayType>>;
  using Value = typename View::non_const_value_type;
  constexpr auto Rank = View::rank();

  if (array.datatype() == DataType::kind<Value>() && array.rank() == Rank) {
    return makeView(array, Value{}, std::integral_constant<int, Rank>{});
  }

  if constexpr (TypeIndex < std::variant_size_v<VariantType<ArrayType>> - 1) {
    return executeMakeView<TypeIndex + 1>(array, makeView);
  } else {
    throw_Exception("ArrayView<" + array.datatype().str() + ", " +
                        std::to_string(array.rank()) +
                        "> is not an alternative in ArrayViewVariant.",
                    Here());
  }
}

template <typename ArrayType>
VariantType<ArrayType> makeViewVariantImpl(ArrayType& array) {
  const auto makeView = [](auto& array, auto value, auto rank) {
    return make_view<decltype(value), decltype(rank)::value>(array);
  };
  return executeMakeView<>(array, makeView);
}

template <typename ArrayType>
VariantType<ArrayType> makeHostViewVariantImpl(ArrayType& array) {
  const auto makeView = [](auto& array, auto value, auto rank) {
    return make_host_view<decltype(value), decltype(rank)::value>(array);
  };
  return executeMakeView<>(array, makeView);
}

template <typename ArrayType>
VariantType<ArrayType> makeDeviceViewVariantImpl(ArrayType& array) {
  const auto makeView = [](auto& array, auto value, auto rank) {
    return make_device_view<decltype(value), decltype(rank)::value>(array);
  };
  return executeMakeView<>(array, makeView);
}

}  // namespace

ArrayViewVariant make_view_variant(Array& array) {
  return makeViewVariantImpl(array);
}

ConstArrayViewVariant make_view_variant(const Array& array) {
  return makeViewVariantImpl(array);
}

ArrayViewVariant make_host_view_variant(Array& array) {
  return makeHostViewVariantImpl(array);
}

ConstArrayViewVariant make_host_view_variant(const Array& array) {
  return makeHostViewVariantImpl(array);
}

ArrayViewVariant make_device_view_variant(Array& array) {
  return makeDeviceViewVariantImpl(array);
}

ConstArrayViewVariant make_device_view_variant(const Array& array) {
  return makeDeviceViewVariantImpl(array);
}

}  // namespace array
}  // namespace atlas
