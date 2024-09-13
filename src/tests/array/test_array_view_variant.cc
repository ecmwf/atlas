/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <type_traits>
#include <variant>

#include "atlas/array.h"
#include "atlas/array/ArrayViewVariant.h"
#include "eckit/utils/Overloaded.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using namespace array;

CASE("test variant assignment") {
  auto array1 = array::ArrayT<float>(2);
  auto array2 = array::ArrayT<double>(2, 3);
  auto array3 = array::ArrayT<int>(2, 3, 4);
  const auto& arrayRef = array1;

  array1.allocateDevice();
  array2.allocateDevice();
  array3.allocateDevice();

  auto view1 = make_view_variant(array1);
  auto view2 = make_view_variant(array2);
  auto view3 = make_view_variant(array3);
  auto view4 = make_view_variant(arrayRef);

  auto hostView1 = make_host_view_variant(array1);
  auto hostView2 = make_host_view_variant(array2);
  auto hostView3 = make_host_view_variant(array3);
  auto hostView4 = make_host_view_variant(arrayRef);

  auto deviceView1 = make_device_view_variant(array1);
  auto deviceView2 = make_device_view_variant(array2);
  auto deviceView3 = make_device_view_variant(array3);
  auto deviceView4 = make_device_view_variant(arrayRef);

  const auto visitVariants = [](auto& var1, auto& var2, auto var3, auto var4) {
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

    std::visit(
        [](auto&& view) {
          using View = std::remove_reference_t<decltype(view)>;
          EXPECT_EQ(View::rank(), 1);
          EXPECT((std::is_same_v<typename View::value_type, const float>));
        },
        var4);
  };

  visitVariants(view1, view2, view3, view4);
  visitVariants(hostView1, hostView2, hostView3, hostView4);
  visitVariants(deviceView1, deviceView2, deviceView3, deviceView4);
}

template <typename View>
constexpr auto Rank() {
  return std::remove_reference_t<View>::rank();
}

CASE("test std::visit") {
  auto array1 = ArrayT<int>(10);
  make_view<int, 1>(array1).assign({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});

  auto array2 = ArrayT<int>(5, 2);
  make_view<int, 2>(array2).assign({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});

  const auto testRank1 = [](auto&& view) {
    using View = std::remove_reference_t<decltype(view)>;
    EXPECT_EQ(View::rank(), 1);
    using Value = typename View::value_type;
    EXPECT((std::is_same_v<Value, int>));

    for (auto i = size_t{0}; i < view.size(); ++i) {
      EXPECT_EQ(view(i), static_cast<Value>(i));
    }
  };

  const auto testRank2 = [](auto&& view) {
    using View = std::remove_reference_t<decltype(view)>;
    EXPECT_EQ(View::rank(), 2);
    using Value = typename View::value_type;
    EXPECT((std::is_same_v<Value, int>));

    auto testValue = int{0};
    for (auto i = idx_t{0}; i < view.shape(0); ++i) {
      for (auto j = idx_t{0}; j < view.shape(1); ++j) {
        EXPECT_EQ(view(i, j), static_cast<Value>(testValue++));
      }
    }
  };

  const auto var1 = make_view_variant(array1);
  const auto var2 = make_view_variant(array2);

  SECTION("demonstrate 'if constexpr' pattern") {
    auto rank1Tested = false;
    auto rank2Tested = false;

    const auto visitor = [&](auto&& view) {
      if constexpr (Rank<decltype(view)>() == 1) {
        testRank1(view);
        rank1Tested = true;
      }
      if constexpr (Rank<decltype(view)>() == 2) {
        testRank2(view);
        rank2Tested = true;
      }
    };

    std::visit(visitor, var1);
    EXPECT(rank1Tested);
    std::visit(visitor, var2);
    EXPECT(rank2Tested);
  }

  SECTION("demonstrate 'overloaded' pattern") {
    // Note, SFINAE can be eliminated using concepts and explicit lambda
    // templates in C++20.
    auto rank1Tested = false;
    auto rank2Tested = false;
    const auto visitor = eckit::Overloaded{
        [&](auto&& view) -> std::enable_if_t<Rank<decltype(view)>() == 1> {
          testRank1(view);
          rank1Tested = true;
        },
        [&](auto&& view) -> std::enable_if_t<Rank<decltype(view)>() == 2> {
          testRank2(view);
          rank2Tested = true;
        },
        [](auto&& view) -> std::enable_if_t<(Rank<decltype(view)>() > 2)> {
          // do nothing.
        }};

    std::visit(visitor, var1);
    EXPECT(rank1Tested);
    std::visit(visitor, var2);
    EXPECT(rank2Tested);
  }
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
