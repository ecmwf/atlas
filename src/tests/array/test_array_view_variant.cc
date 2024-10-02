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

  const auto hostView1 = make_host_view_variant(array1);
  const auto hostView2 = make_host_view_variant(array2);
  const auto hostView3 = make_host_view_variant(array3);
  const auto hostView4 = make_host_view_variant(arrayRef);

  auto deviceView1 = make_device_view_variant(array1);
  auto deviceView2 = make_device_view_variant(array2);
  auto deviceView3 = make_device_view_variant(array3);
  auto deviceView4 = make_device_view_variant(arrayRef);

  auto view = make_view<float, 1>(array1);

  const auto visitVariants = [](auto& var1, auto& var2, auto var3, auto var4) {
    std::visit(
        [](auto&& view) {
          EXPECT((is_rank<1>(view)));
          EXPECT((is_value_type<float>(view)));
          EXPECT((is_non_const_value_type<float>(view)));
        },
        var1);

    std::visit(
        [](auto&& view) {
          EXPECT((is_rank<2>(view)));
          EXPECT((is_value_type<double>(view)));
          EXPECT((is_non_const_value_type<double>(view)));
        },
        var2);

    std::visit(
        [](auto&& view) {
          EXPECT((is_rank<3>(view)));
          EXPECT((is_value_type<int>(view)));
          EXPECT((is_non_const_value_type<int>(view)));
        },
        var3);

    std::visit(
        [](auto&& view) {
          EXPECT((is_rank<1>(view)));
          EXPECT((is_value_type<const float>(view)));
          EXPECT((is_non_const_value_type<float>(view)));
        },
        var4);
  };

  visitVariants(view1, view2, view3, view4);
  visitVariants(hostView1, hostView2, hostView3, hostView4);
  visitVariants(deviceView1, deviceView2, deviceView3, deviceView4);
}

CASE("test std::visit") {
  auto array1 = ArrayT<int>(10);
  make_view<int, 1>(array1).assign({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});

  auto array2 = ArrayT<int>(5, 2);
  make_view<int, 2>(array2).assign({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});

  const auto var1 = make_view_variant(array1);
  const auto var2 = make_view_variant(array2);
  auto rank1Tested = false;
  auto rank2Tested = false;

  const auto visitor = [&](auto&& view) {
    if constexpr (is_rank<1>(view)) {
      EXPECT((is_value_type<int>(view)));
      auto testValue = int{0};
      for (auto i = size_t{0}; i < view.size(); ++i) {
        const auto value = view(i);
        EXPECT_EQ(value, static_cast<decltype(value)>(testValue++));
      }
      rank1Tested = true;
    } else if constexpr (is_rank<2>(view)) {
      EXPECT((is_value_type<int>(view)));
      auto testValue = int{0};
      for (auto i = idx_t{0}; i < view.shape(0); ++i) {
        for (auto j = idx_t{0}; j < view.shape(1); ++j) {
          const auto value = view(i, j);
          EXPECT_EQ(value, static_cast<decltype(value)>(testValue++));
        }
      }
      rank2Tested = true;
    } else {
      // Test should not reach here.
      EXPECT(false);
    }
  };

  std::visit(visitor, var1);
  EXPECT(rank1Tested);
  std::visit(visitor, var2);
  EXPECT(rank2Tested);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
