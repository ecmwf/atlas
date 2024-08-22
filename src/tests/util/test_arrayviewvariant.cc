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
#include "atlas/util/ArrayViewVariant.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using namespace array;
using namespace util;

CASE("test visit") {
  auto arr1 = array::ArrayT<float>(2);
  auto arr2 = array::ArrayT<double>(2, 3);
  auto arr3 = array::ArrayT<int>(2, 3, 4);

  using ValueList = Values<float, double, int>;
  using RankList = Ranks<1, 2, 3>;

  const auto var1 = make_array_view_variant<ValueList, RankList>(arr1);
  const auto var2 = make_array_view_variant<ValueList, RankList>(arr2);
  const auto var3 = make_array_view_variant<ValueList, RankList>(arr3);

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

CASE("test array view data") {
  auto arr = ArrayT<int>(10);
  make_view<int, 1>(arr).assign({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});

  const auto& arrRef = arr;
  const auto var =
      make_array_view_variant<Values<int, double>, Ranks<1, 2>>(arrRef);

  std::visit(
      [](auto&& view) {
        using View = std::remove_reference_t<decltype(view)>;
        EXPECT_EQ(View::rank(), 1);
        if constexpr (View::rank() == 1) {
          for (auto i = idx_t{0}; i < view.size(); ++i) {
            EXPECT_EQ(view(i), i);
          }
        }
      },
      var);
}

CASE("test instantiation") {
  auto arr = array::ArrayT<double>(1);
  const auto constArr = array::ArrayT<double>(1);

  SECTION("default variants") {
    auto var = make_array_view_variant(arr);
    auto constVar = make_array_view_variant(constArr);

    using VarType = std::variant<ArrayView<float, 1>, ArrayView<float, 2>,
                                 ArrayView<float, 3>, ArrayView<double, 1>,
                                 ArrayView<double, 2>, ArrayView<double, 3>>;

    using ConstVarType =
        std::variant<ArrayView<const float, 1>, ArrayView<const float, 2>,
                     ArrayView<const float, 3>, ArrayView<const double, 1>,
                     ArrayView<const double, 2>, ArrayView<const double, 3>>;

    EXPECT((std::is_same_v<decltype(var), VarType>));
    EXPECT((std::is_same_v<decltype(constVar), ConstVarType>));
  }

  SECTION("customised variants") {
    using ValueList = Values<int, double>;
    using RankList = Ranks<1>;

    auto var = make_array_view_variant<ValueList, RankList>(arr);
    auto constVar = make_array_view_variant<ValueList, RankList>(constArr);

    using VarType = std::variant<ArrayView<int, 1>, ArrayView<double, 1>>;

    using ConstVarType =
        std::variant<ArrayView<const int, 1>, ArrayView<const double, 1>>;

    EXPECT((std::is_same_v<decltype(var), VarType>));
    EXPECT((std::is_same_v<decltype(constVar), ConstVarType>));
  }

  SECTION("mismatched array and variant") {
    // Array is of value-type double.
    using ValueList = Values<int, float>;
    using RankList = Ranks<1>;

    EXPECT_THROWS((make_array_view_variant<ValueList, RankList>(arr)));
    EXPECT_THROWS((make_array_view_variant<ValueList, RankList>(constArr)));
  }
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
