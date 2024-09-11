/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>
#include <type_traits>
#include <variant>

#include "atlas/array.h"
#include "atlas/array/ArrayViewVariant.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using namespace array;

CASE("test visit") {
  auto arr1 = array::ArrayT<float>(2);
  auto arr2 = array::ArrayT<double>(2, 3);
  auto arr3 = array::ArrayT<int>(2, 3, 4);

  const auto var1 = make_view_variant(arr1);
  const auto var2 = make_view_variant(arr2);
  const auto var3 = make_view_variant(arr3);

  std::visit(
      [](auto&& view) {
        using View = std::remove_reference_t<decltype(view)>;
        EXPECT_EQ(View::rank(), 1);
        EXPECT((std::is_same_v<typename View::value_type, float>));
      },
      var1);

  std::visit(
      [](auto&& view) {
        using View = std::remove_reference_t<decltype(view)>;
        EXPECT_EQ(View::rank(), 2);
        EXPECT((std::is_same_v<typename View::value_type, double>));
      },
      var2);

  std::visit(
      [](auto&& view) {
        using View = std::remove_reference_t<decltype(view)>;
        EXPECT_EQ(View::rank(), 3);
        EXPECT((std::is_same_v<typename View::value_type, int>));
      },
      var3);
}

template <typename View>
constexpr auto Rank = std::decay_t<View>::rank();

CASE("test array view data") {
  auto arr = ArrayT<int>(10);
  make_view<int, 1>(arr).assign({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});

  const auto& arrRef = arr;
  const auto var = make_view_variant(arrRef);

  const auto visitor = Overloaded{
      [](auto&& view) -> std::enable_if_t<Rank<decltype(view)> == 1> {
        using View = std::decay_t<decltype(view)>;
        EXPECT_EQ(View::rank(), 1);
        using Value = typename View::value_type;
        EXPECT((std::is_same_v<Value, const int>));

        for (auto i = size_t{0}; i < view.size(); ++i) {
          EXPECT_EQ(view(i), static_cast<Value>(i));
        }
      },
      [](auto&& view) -> std::enable_if_t<Rank<decltype(view)> != 1> {
        // do nothing.
      }};

  std::visit(visitor, var);
}

CASE("test instantiation") {
  auto arr = array::ArrayT<double>(1);
  const auto constArr = array::ArrayT<double>(1);

  //  SECTION("default variants") {
  //    auto var = make_view_variant(arr);
  //    auto constVar = make_view_variant(constArr);

  //    using VarType = std::variant<ArrayView<float, 1>, ArrayView<float, 2>,
  //                                 ArrayView<float, 3>, ArrayView<double, 1>,
  //                                 ArrayView<double, 2>, ArrayView<double,
  //                                 3>>;

  //    using ConstVarType =
  //        std::variant<ArrayView<const float, 1>, ArrayView<const float, 2>,
  //                     ArrayView<const float, 3>, ArrayView<const double, 1>,
  //                     ArrayView<const double, 2>, ArrayView<const double,
  //                     3>>;

  //    EXPECT((std::is_same_v<decltype(var), VarType>));
  //    EXPECT((std::is_same_v<decltype(constVar), ConstVarType>));
  //  }
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
