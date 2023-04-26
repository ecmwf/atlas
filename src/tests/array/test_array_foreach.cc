/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <chrono>
#include <type_traits>

#include "atlas/array.h"
#include "atlas/array/MakeView.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/array/helpers/ArraySlicer.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::array;
using namespace atlas::array::helpers;

namespace atlas {
namespace test {

CASE("test_array_foreach_1_view") {

  const auto arr = ArrayT<double>(2, 3);
  const auto view = make_view<double, 2>(arr);

  // Test slice shapes.

  const auto loopFunctorDim0 = [](auto& slice) {
    EXPECT_EQUAL(slice.rank(), 1);
    EXPECT_EQUAL(slice.shape(0), 3);
  };
  ArrayForEach<0>::apply(std::make_tuple(view), loopFunctorDim0);

  const auto loopFunctorDim1 = [](auto& slice) {
    EXPECT_EQUAL(slice.rank(), 1);
    EXPECT_EQUAL(slice.shape(0), 2);
  };
  ArrayForEach<1>::apply(std::make_tuple(view), loopFunctorDim1);

  // Test that slice resolves to double.

  const auto loopFunctorDimAll = [](auto& slice) {
    static_assert(std::is_convertible_v<decltype(slice), const double&>);
  };
  ArrayForEach<0, 1>::apply(std::make_tuple(view), loopFunctorDimAll);

  // Test ghost functionality.

  auto ghost = ArrayT<int>(2);
  auto ghostView = make_view<int, 1>(ghost);
  ghostView.assign({0, 1});

  auto count = int {};
  const auto countNonGhosts = [&count](auto&...) { ++count; };
  ArrayForEach<0>::apply(std::make_tuple(view), countNonGhosts, ghostView,
                         array::helpers::detail::sequencedConf());
  EXPECT_EQ(count, 1);

  count = 0;
  const auto ghostWrap = [&ghostView](idx_t idx, auto&...) {
    // Wrap ghostView to use correct number of indices.
    return ghostView(idx);
  };
  ArrayForEach<0, 1>::apply(std::make_tuple(view), countNonGhosts, ghostWrap,
                            array::helpers::detail::sequencedConf());
  EXPECT_EQ(count, 3);
}

CASE("test_array_foreach_2_views") {

  const auto arr1 = ArrayT<double>(2, 3);
  const auto view1 = make_view<double, 2>(arr1);

  const auto arr2 = ArrayT<double>(2, 3, 4);
  const auto view2 = make_view<double, 3>(arr2);

  // Test slice shapes.

  const auto loopFunctorDim0 = [](auto& slice1, auto& slice2) {
    EXPECT_EQUAL(slice1.rank(), 1);
    EXPECT_EQUAL(slice1.shape(0), 3);

    EXPECT_EQUAL(slice2.rank(), 2);
    EXPECT_EQUAL(slice2.shape(0), 3);
    EXPECT_EQUAL(slice2.shape(1), 4);
  };
  ArrayForEach<0>::apply(std::make_tuple(view1, view2), loopFunctorDim0);

  const auto loopFunctorDim1 = [](auto& slice1, auto& slice2) {
    EXPECT_EQUAL(slice1.rank(), 1);
    EXPECT_EQUAL(slice1.shape(0), 2);

    EXPECT_EQUAL(slice2.rank(), 2);
    EXPECT_EQUAL(slice2.shape(0), 2);
    EXPECT_EQUAL(slice2.shape(1), 4);
  };
  ArrayForEach<1>::apply(std::make_tuple(view1, view2), loopFunctorDim1);

  // Test that slice resolves to double.

  const auto loopFunctorDimAll = [](auto& slice2) {
    static_assert(std::is_convertible_v<decltype(slice2), const double&>);
  };
  ArrayForEach<0, 1, 2>::apply(std::make_tuple(view2), loopFunctorDimAll);

  // Test ghost functionality.

  auto ghost = ArrayT<int>(2);
  auto ghostView = make_view<int, 1>(ghost);
  ghostView.assign({0, 1});

  auto count = int {};
  const auto countNonGhosts = [&count](auto&...) { ++count; };
  ArrayForEach<0>::apply(std::make_tuple(view2), countNonGhosts, ghostView,
                         array::helpers::detail::sequencedConf());
  EXPECT_EQ(count, 1);

  count = 0;
  const auto ghostWrap = [&ghostView](idx_t idx, auto&...) {
    // Wrap ghostView to use correct number of indices.
    return ghostView(idx);
  };
  ArrayForEach<0, 1>::apply(std::make_tuple(view2), countNonGhosts, ghostWrap,
                            array::helpers::detail::sequencedConf());
  EXPECT_EQ(count, 3);

  count = 0;
  ArrayForEach<0, 1, 2>::apply(std::make_tuple(view2), countNonGhosts,
                               ghostWrap,
                               array::helpers::detail::sequencedConf());
  EXPECT_EQ(count, 12);
}

CASE("test_array_foreach_3_views") {

  const auto arr1 = ArrayT<double>(2, 3);
  const auto view1 = make_view<double, 2>(arr1);

  const auto arr2 = ArrayT<double>(2, 3, 4);
  const auto view2 = make_view<double, 3>(arr2);

  const auto arr3 = ArrayT<double>(2, 3, 4, 5);
  const auto view3 = make_view<double, 4>(arr3);

  // Test slice shapes.

  const auto loopFunctorDim0 = [](auto& slice1, auto& slice2, auto& slice3) {
    EXPECT_EQUAL(slice1.rank(), 1);
    EXPECT_EQUAL(slice1.shape(0), 3);

    EXPECT_EQUAL(slice2.rank(), 2);
    EXPECT_EQUAL(slice2.shape(0), 3);
    EXPECT_EQUAL(slice2.shape(1), 4);

    EXPECT_EQUAL(slice3.rank(), 3);
    EXPECT_EQUAL(slice3.shape(0), 3);
    EXPECT_EQUAL(slice3.shape(1), 4);
    EXPECT_EQUAL(slice3.shape(2), 5);
  };
  ArrayForEach<0>::apply(std::make_tuple(view1, view2, view3), loopFunctorDim0);

  const auto loopFunctorDim1 = [](auto& slice1, auto& slice2, auto& slice3) {
    EXPECT_EQUAL(slice1.rank(), 1);
    EXPECT_EQUAL(slice1.shape(0), 2);

    EXPECT_EQUAL(slice2.rank(), 2);
    EXPECT_EQUAL(slice2.shape(0), 2);
    EXPECT_EQUAL(slice2.shape(1), 4);

    EXPECT_EQUAL(slice3.rank(), 3);
    EXPECT_EQUAL(slice3.shape(0), 2);
    EXPECT_EQUAL(slice3.shape(1), 4);
    EXPECT_EQUAL(slice3.shape(2), 5);
  };
  ArrayForEach<1>::apply(std::make_tuple(view1, view2, view3), loopFunctorDim1);

  // Test that slice resolves to double.

  const auto loopFunctorDimAll = [](auto& slice3) {
    static_assert(std::is_convertible_v<decltype(slice3), const double&>);
  };
  ArrayForEach<0, 1, 2, 3>::apply(std::make_tuple(view3), loopFunctorDimAll);

  // Test ghost functionality.

  auto ghost = ArrayT<int>(2);
  auto ghostView = make_view<int, 1>(ghost);
  ghostView.assign({0, 1});

  auto count = int {};
  const auto countNonGhosts = [&count](auto&...) { ++count; };
  ArrayForEach<0>::apply(std::make_tuple(view3), countNonGhosts, ghostView,
                         array::helpers::detail::sequencedConf());
  EXPECT_EQ(count, 1);

  count = 0;
  const auto ghostWrap = [&ghostView](idx_t idx, auto&...) {
    // Wrap ghostView to use correct number of indices.
    return ghostView(idx);
  };
  ArrayForEach<0, 1>::apply(std::make_tuple(view3), countNonGhosts, ghostWrap,
                            array::helpers::detail::sequencedConf());
  EXPECT_EQ(count, 3);

  count = 0;
  ArrayForEach<0, 1, 2>::apply(std::make_tuple(view3), countNonGhosts,
                               ghostWrap,
                               array::helpers::detail::sequencedConf());
  EXPECT_EQ(count, 12);

  count = 0;
  ArrayForEach<0, 1, 2, 3>::apply(std::make_tuple(view3), countNonGhosts,
                                  ghostWrap,
                                  array::helpers::detail::sequencedConf());
  EXPECT_EQ(count, 60);
}

CASE("test_array_foreach_data_integrity") {

  auto arr1 = ArrayT<double>(200, 3);
  auto view1 = make_view<double, 2>(arr1);

  auto arr2 = ArrayT<double>(200, 3, 4);
  auto view2 = make_view<double, 3>(arr2);

  for (auto idx = size_t{}; idx < arr1.size(); ++idx) {
    static_cast<double*>(arr1.data())[idx] = idx;
  }

  for (auto idx = size_t{}; idx < arr2.size(); ++idx) {
    static_cast<double*>(arr2.data())[idx] = idx;
  }

  const auto scaleDataDim0 = [](auto& slice1, auto& slice2) {

    static_assert(std::is_convertible_v<decltype(slice1), double&>);
    slice1 *= 2.;

    const auto scaleDataDim1 = [](auto& slice) {

      static_assert(std::is_convertible_v<decltype(slice), double&>);
      slice *= 3.;
    };
    ArrayForEach<0>::apply(std::make_tuple(slice2), scaleDataDim1,
                           array::helpers::detail::sequencedConf());
  };
  ArrayForEach<0, 1>::apply(std::make_tuple(view1, view2), scaleDataDim0);

  for (auto idx = size_t{}; idx < arr1.size(); ++idx) {
    EXPECT_EQ(static_cast<double*>(arr1.data())[idx], 2. * idx);
  }

  for (auto idx = size_t{}; idx < arr2.size(); ++idx) {
    EXPECT_EQ(static_cast<double*>(arr2.data())[idx], 3. * idx);
  }
}

template <typename IterationMethod, typename Operation>
void timeLoop(const IterationMethod& iterationMethod,
              const Operation& operation, const std::string& output) {
  const auto start = std::chrono::steady_clock::now();
  iterationMethod(operation);
  const auto end = std::chrono::steady_clock::now();
  const auto duration = std::chrono::duration<double>{end - start};

  Log::info() << "Elapsed time: " + output + "= " << duration.count() << "s\n";
}

CASE("test_array_foreach_performance") {

  auto arr1 = ArrayT<double>(50000, 100);
  auto view1 = make_view<double, 2>(arr1);

  auto arr2 = ArrayT<double>(50000, 100);
  auto view2 = make_view<double, 2>(arr2);

  for (auto idx = size_t{}; idx < arr2.size(); ++idx) {
    static_cast<double*>(arr2.data())[idx] = 2 * idx + 1;
  }

  auto arr3 = ArrayT<double>(50000, 100);
  auto view3 = make_view<double, 2>(arr2);

  for (auto idx = size_t{}; idx < arr3.size(); ++idx) {
    static_cast<double*>(arr3.data())[idx] = 3 * idx + 1;
  }

  const auto add = [](double& a1, const double& a2,
                      const double& a3) { a1 = a2 + a3; };

  const auto trig = [](double& a1, const double& a2,
                       const double& a3) { a1 = std::sin(a2) + std::cos(a3); };

  const auto rawPointer = [&](const auto& operation) {
    auto* p1 = static_cast<double*>(arr1.data());
    const auto* p2 = static_cast<double*>(arr2.data());
    const auto* p3 = static_cast<double*>(arr3.data());
    for (auto idx = size_t{}; idx < arr1.size(); ++idx) {
      operation(p1[idx], p2[idx], p3[idx]);
    }
  };

  const auto ijLoop = [&](const auto& operation) {
    for (auto i = idx_t{}; i < view1.shape(0); ++i) {
      for (auto j = idx_t{}; j < view1.shape(1); ++j) {
        operation(view1(i, j), view2(i, j), view3(i, j));
      }
    }
  };

  const auto jiLoop = [&](const auto& operation) {
    for (auto j = idx_t{}; j < view1.shape(1); ++j) {
      for (auto i = idx_t{}; i < view1.shape(0); ++i) {
        operation(view1(i, j), view2(i, j), view3(i, j));
      }
    }
  };

  const auto forEachCol = [&](const auto& operation) {
    const auto function = [&](auto& slice1, auto& slice2, auto& slice3) {
      for (auto idx = idx_t{}; idx < slice1.shape(0); ++idx) {
        operation(slice1(idx), slice2(idx), slice3(idx));
      }
    };
    ArrayForEach<0>::apply(std::make_tuple(view1, view2, view3), function,
                           array::helpers::detail::sequencedConf());
  };

  const auto forEachLevel = [&](const auto& operation) {
    const auto function = [&](auto& slice1, auto& slice2, auto& slice3) {
      for (auto idx = idx_t{}; idx < slice1.shape(0); ++idx) {
        operation(slice1(idx), slice2(idx), slice3(idx));
      }
    };
    ArrayForEach<1>::apply(std::make_tuple(view1, view2, view3), function,
                           array::helpers::detail::sequencedConf());
  };

  const auto forEachAll = [&](const auto& operation) {
    ArrayForEach<0, 1>::apply(std::make_tuple(view1, view2, view3), operation,
                              array::helpers::detail::sequencedConf());
  };

  timeLoop(rawPointer, add, "Addition; raw pointer             ");
  timeLoop(ijLoop, add, "Addition; for loop (i, j)         ");
  timeLoop(jiLoop, add, "Addition; for loop (j, i)         ");
  timeLoop(forEachCol, add, "Addition; for each (columns)      ");
  timeLoop(forEachLevel, add, "Addition; for each (levels)       ");
  timeLoop(forEachAll, add, "Addition; for each (all elements) ");
  timeLoop(rawPointer, trig, "Trig    ; raw pointer             ");
  timeLoop(ijLoop, trig, "Trig    ; for loop (i, j)         ");
  timeLoop(jiLoop, trig, "Trig    ; for loop (j, i)         ");
  timeLoop(forEachCol, trig, "Trig    ; for each (columns)      ");
  timeLoop(forEachLevel, trig, "Trig    ; for each (levels)       ");
  timeLoop(forEachAll, trig, "Trig    ; for each (all elements) ");
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
