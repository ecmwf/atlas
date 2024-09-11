/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/ArrayViewVariant.h"

#include <string>

#include "atlas/runtime/Exception.h"

namespace atlas {
namespace array {

using namespace detail;

namespace {

// Match array.rank() and array.datatype() to variant types. Return result of
// makeView on a successfull pattern match.
template <size_t TypeIndex>
struct MakeViewHelper {
  template <typename ArrayType, typename MakeView>
  static ArrayViewVariant apply(ArrayType& array, const MakeView& makeView) {
    using View = std::variant_alternative_t<TypeIndex, ArrayViewVariant>;
    using Value = typename View::value_type;
    constexpr auto Rank = View::rank();

    if (array.datatype() == DataType::kind<Value>() && array.rank() == Rank) {
      return makeView(array, Value{}, std::integral_constant<idx_t, Rank>{});
    }
    return MakeViewHelper<TypeIndex + 1>::apply(array, makeView);
  }
};

// End recursion.
template <>
struct MakeViewHelper<std::variant_size_v<ArrayViewVariant>> {
  template <typename ArrayType, typename MakeView>
  static ArrayViewVariant apply(ArrayType& array, const MakeView& makeView) {
    throw_Exception("ArrayView<" + array.datatype().str() + ", " +
                        std::to_string(array.rank()) +
                        "> is not an alternative in ArrayViewVariant.",
                    Here());
  }
};

template <typename ArrayType>
ArrayViewVariant makeViewVariantImpl(ArrayType& array) {
  const auto makeView = [](auto& array, auto value, auto rank) {
    return make_view<decltype(value), decltype(rank)::value>(array);
  };
  return MakeViewHelper<0>::apply(array, makeView);
}

template <typename ArrayType>
ArrayViewVariant makeHostViewVariantImpl(ArrayType& array) {
  const auto makeView = [](auto& array, auto value, auto rank) {
    return make_host_view<decltype(value), decltype(rank)::value>(array);
  };
  return MakeViewHelper<0>::apply(array, makeView);
}

template <typename ArrayType>
ArrayViewVariant makeDeviceViewVariantImpl(ArrayType& array) {
  const auto makeView = [](auto& array, auto value, auto rank) {
    return make_device_view<decltype(value), decltype(rank)::value>(array);
  };
  return MakeViewHelper<0>::apply(array, makeView);
}

}  // namespace

ArrayViewVariant make_view_variant(Array& array) {
  return makeViewVariantImpl(array);
}

ArrayViewVariant make_view_variant(const Array& array) {
  return makeViewVariantImpl(array);
}

ArrayViewVariant make_host_view_variant(Array& array) {
  return makeHostViewVariantImpl(array);
}

ArrayViewVariant make_host_view_variant(const Array& array) {
  return makeHostViewVariantImpl(array);
}

ArrayViewVariant make_devive_view_variant(Array& array) {
  return makeDeviceViewVariantImpl(array);
}

ArrayViewVariant make_device_view_variant(const Array& array) {
  return makeDeviceViewVariantImpl(array);
}

}  // namespace array
}  // namespace atlas
