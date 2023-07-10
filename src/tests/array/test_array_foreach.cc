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

  const auto loopFunctorDim0 = [](auto&& slice) {
    EXPECT_EQ(slice.rank(), 1);
    EXPECT_EQ(slice.shape(0), 3);
  };
  ArrayForEach<0>::apply(std::tie(view), loopFunctorDim0);

  const auto loopFunctorDim1 = [](auto&& slice) {
    EXPECT_EQ(slice.rank(), 1);
    EXPECT_EQ(slice.shape(0), 2);
  };
  ArrayForEach<1>::apply(std::tie(view), loopFunctorDim1);

  // Test that slice resolves to double.

  const auto loopFunctorDimAll = [](auto&& slice) {
    static_assert(std::is_convertible_v<decltype(slice), const double&>);
  };
  ArrayForEach<0, 1>::apply(std::tie(view), loopFunctorDimAll);

  // Test ghost functionality.

  auto ghost = ArrayT<int>(2);
  auto ghostView = make_view<int, 1>(ghost);
  ghostView.assign({0, 1});

  auto count = int {};
  const auto countNonGhosts = [&count](auto&&...) { ++count; };
  ArrayForEach<0>::apply(execution::seq, std::tie(view), ghostView, countNonGhosts);
  EXPECT_EQ(count, 1);

  count = 0;
  const auto ghostWrap = [&ghostView](idx_t idx, auto&&...) {
    // Wrap ghostView to use correct number of indices.
    return ghostView(idx);
  };
  ArrayForEach<0, 1>::apply(execution::seq, std::tie(view), ghostWrap, countNonGhosts);
  EXPECT_EQ(count, 3);
}

CASE("test_array_foreach_2_views") {

  const auto arr1 = ArrayT<double>(2, 3);
  const auto view1 = make_view<double, 2>(arr1);

  const auto arr2 = ArrayT<double>(2, 3, 4);
  const auto view2 = make_view<double, 3>(arr2);

  // Test slice shapes.

  const auto loopFunctorDim0 = [](auto&& slice1, auto&& slice2) {
    EXPECT_EQ(slice1.rank(), 1);
    EXPECT_EQ(slice1.shape(0), 3);

    EXPECT_EQ(slice2.rank(), 2);
    EXPECT_EQ(slice2.shape(0), 3);
    EXPECT_EQ(slice2.shape(1), 4);
  };
  ArrayForEach<0>::apply(std::tie(view1, view2), loopFunctorDim0);

  const auto loopFunctorDim1 = [](auto&& slice1, auto&& slice2) {
    EXPECT_EQ(slice1.rank(), 1);
    EXPECT_EQ(slice1.shape(0), 2);

    EXPECT_EQ(slice2.rank(), 2);
    EXPECT_EQ(slice2.shape(0), 2);
    EXPECT_EQ(slice2.shape(1), 4);
  };
  ArrayForEach<1>::apply(std::tie(view1, view2), loopFunctorDim1);

  // Test that slice resolves to double.

  const auto loopFunctorDimAll = [](auto&& slice2) {
    static_assert(std::is_convertible_v<decltype(slice2), const double&>);
  };
  ArrayForEach<0, 1, 2>::apply(std::tie(view2), loopFunctorDimAll);

  // Test ghost functionality.

  auto ghost = ArrayT<int>(2);
  auto ghostView = make_view<int, 1>(ghost);
  ghostView.assign({0, 1});

  auto count = int {};
  const auto countNonGhosts = [&count](auto&&...) { ++count; };
  ArrayForEach<0>::apply(execution::seq, std::tie(view2), ghostView, countNonGhosts);
  EXPECT_EQ(count, 1);

  count = 0;
  const auto ghostWrap = [&ghostView](idx_t idx, auto&&...) {
    // Wrap ghostView to use correct number of indices.
    return ghostView(idx);
  };
  ArrayForEach<0, 1>::apply(execution::seq, std::tie(view2), ghostWrap, countNonGhosts);
  EXPECT_EQ(count, 3);

  count = 0;
  ArrayForEach<0, 1, 2>::apply(execution::seq, std::tie(view2), ghostWrap, countNonGhosts);
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

  const auto loopFunctorDim0 = [](auto&& slice1, auto&& slice2, auto&& slice3) {
    EXPECT_EQ(slice1.rank(), 1);
    EXPECT_EQ(slice1.shape(0), 3);

    EXPECT_EQ(slice2.rank(), 2);
    EXPECT_EQ(slice2.shape(0), 3);
    EXPECT_EQ(slice2.shape(1), 4);

    EXPECT_EQ(slice3.rank(), 3);
    EXPECT_EQ(slice3.shape(0), 3);
    EXPECT_EQ(slice3.shape(1), 4);
    EXPECT_EQ(slice3.shape(2), 5);
  };
  ArrayForEach<0>::apply(std::tie(view1, view2, view3), loopFunctorDim0);

  const auto loopFunctorDim1 = [](auto&& slice1, auto&& slice2, auto&& slice3) {
    EXPECT_EQ(slice1.rank(), 1);
    EXPECT_EQ(slice1.shape(0), 2);

    EXPECT_EQ(slice2.rank(), 2);
    EXPECT_EQ(slice2.shape(0), 2);
    EXPECT_EQ(slice2.shape(1), 4);

    EXPECT_EQ(slice3.rank(), 3);
    EXPECT_EQ(slice3.shape(0), 2);
    EXPECT_EQ(slice3.shape(1), 4);
    EXPECT_EQ(slice3.shape(2), 5);
  };
  ArrayForEach<1>::apply(std::tie(view1, view2, view3), loopFunctorDim1);

  // Test that slice resolves to double.

  const auto loopFunctorDimAll = [](auto&& slice3) {
    static_assert(std::is_convertible_v<decltype(slice3), const double&>);
  };
  ArrayForEach<0, 1, 2, 3>::apply(std::tie(view3), loopFunctorDimAll);

  // Test ghost functionality.

  auto ghost = ArrayT<int>(2);
  auto ghostView = make_view<int, 1>(ghost);
  ghostView.assign({0, 1});

  auto count = int {};
  const auto countNonGhosts = [&count](auto&&...) { ++count; };
  ArrayForEach<0>::apply(execution::seq, std::tie(view3), ghostView, countNonGhosts);
  EXPECT_EQ(count, 1);

  count = 0;
  const auto ghostWrap = [&ghostView](idx_t idx, auto&&...) {
    // Wrap ghostView to use correct number of indices.
    return ghostView(idx);
  };
  ArrayForEach<0, 1>::apply(execution::seq, std::tie(view3), ghostWrap, countNonGhosts);
  EXPECT_EQ(count, 3);

  count = 0;
  ArrayForEach<0, 1, 2>::apply(execution::seq, std::tie(view3), ghostWrap, countNonGhosts);
  EXPECT_EQ(count, 12);

  count = 0;
  ArrayForEach<0, 1, 2, 3>::apply(execution::seq, std::tie(view3), ghostWrap, countNonGhosts);
  EXPECT_EQ(count, 60);
}


CASE("test_array_foreach_forwarding") {

  const auto arr1 = ArrayT<double>(2, 3);
  const auto view1 = make_view<double, 2>(arr1);

  auto arr2 = ArrayT<double>(2, 3, 4);
  auto view2 = make_view<double, 3>(arr2);

  const auto loopFunctorDim0 = [](auto&& slice1, auto&& slice2) {
    EXPECT_EQUAL(slice1.rank(), 1);
    EXPECT_EQUAL(slice1.shape(0), 3);

    EXPECT_EQUAL(slice2.rank(), 2);
    EXPECT_EQUAL(slice2.shape(0), 3);
    EXPECT_EQUAL(slice2.shape(1), 4);
  };

  ArrayForEach<0>::apply(std::make_tuple(view1, view2), loopFunctorDim0);
  ArrayForEach<0>::apply(std::tie(view1, view2), loopFunctorDim0);
  ArrayForEach<0>::apply(std::forward_as_tuple(view1, view2), loopFunctorDim0);
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

  const auto scaleDataDim0 = [](auto&& slice1, auto&& slice2) {

    static_assert(std::is_convertible_v<decltype(slice1), double&>);
    slice1 *= 2.;

    const auto scaleDataDim1 = [](auto&& slice) {

      static_assert(std::is_convertible_v<decltype(slice), double&>);
      slice *= 3.;
    };
    ArrayForEach<0>::apply(execution::seq, std::tie(slice2), scaleDataDim1);
  };
  ArrayForEach<0, 1>::apply(std::tie(view1, view2), scaleDataDim0);

  for (auto idx = size_t{}; idx < arr1.size(); ++idx) {
    EXPECT_EQ(static_cast<double*>(arr1.data())[idx], 2. * idx);
  }

  for (auto idx = size_t{}; idx < arr2.size(); ++idx) {
    EXPECT_EQ(static_cast<double*>(arr2.data())[idx], 3. * idx);
  }
}

template <typename IterationMethod, typename Operation>
double timeLoop(const IterationMethod& iterationMethod, int num_iter, int num_first,
                const Operation& operation, double baseline, const std::string& output) {
  double time{0};
  for (int i=0; i<num_first + num_iter; ++i) {
    const auto start = std::chrono::steady_clock::now();
    iterationMethod(operation);
    const auto stop = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration<double>{stop - start};
    if (i>=num_first) {
      time += duration.count();
    }
  }
  time /= double(num_iter);
  Log::info() << "Elapsed time: " + output + "= " << time << "s";
  if (baseline != 0) {
    Log::info() << "\t;   relative to baseline : " << 100.*time/baseline << "%";
  }
  Log::info() << std::endl;
  return time;
}

CASE("test_array_foreach_performance") {
  int ni = 50000;
  int nj = 100;
  int num_iter = 20;
  int num_first = 3;

  if( ATLAS_ARRAYVIEW_BOUNDS_CHECKING ) {
    ni = 5000;
    nj = 20;
    num_iter = 1;
    num_first = 0;
    Log::info() << "WARNING: Following timings contain very expensive bounds checking. Compile with -DENABLE_BOUNDSCHECKING=OFF for fair comparison" << std::endl;
  }



  auto arr1 = ArrayT<double>(ni, nj);
  auto view1 = make_view<double, 2>(arr1);

  auto arr2 = ArrayT<double>(ni, nj);
  auto view2 = make_view<double, 2>(arr2);

  for (auto idx = size_t{}; idx < arr2.size(); ++idx) {
    static_cast<double*>(arr2.data())[idx] = 2 * idx + 1;
  }

  auto arr3 = ArrayT<double>(ni, nj);
  auto view3 = make_view<double, 2>(arr2);


  for (auto idx = size_t{}; idx < arr3.size(); ++idx) {
    static_cast<double*>(arr3.data())[idx] = 3 * idx + 1;
  }

  const auto add = [](double& a1, const double& a2,
                      const double& a3) { a1 = a2 + a3; };

  const auto trig = [](double& a1, const double& a2,
                       const double& a3) { a1 = std::sin(a2) + std::cos(a3); };

  const auto rawPointer = [&](const auto& operation) {
    const size_t size = arr1.size();
    auto* p1 = view1.data();
    const auto* p2 = view2.data();
    const auto* p3 = view3.data();
    for (size_t idx = 0; idx < size; ++idx) {
      operation(p1[idx], p2[idx], p3[idx]);
    }
  };

  const auto ijLoop = [&](const auto& operation) {
    const idx_t ni = view1.shape(0);
    const idx_t nj = view1.shape(1);
    for (idx_t i = 0; i < ni; ++i) {
      for (idx_t j = 0; j < nj; ++j) {
        operation(view1(i, j), view2(i, j), view3(i, j));
      }
    }
  };

  const auto jiLoop = [&](const auto& operation) {
    const idx_t ni = view1.shape(0);
    const idx_t nj = view1.shape(1);
    for (idx_t j = 0; j < nj; ++j) {
      for (idx_t i = 0; i < ni; ++i) {
        operation(view1(i, j), view2(i, j), view3(i, j));
      }
    }
  };

  const auto forEachCol = [&](const auto& operation) {
    const auto function = [&](auto&& slice1, auto&& slice2, auto&& slice3) {
      const idx_t size = slice1.shape(0);
      for (idx_t idx = 0; idx < size; ++idx) {
        operation(slice1(idx), slice2(idx), slice3(idx));
      }
    };
    ArrayForEach<0>::apply(execution::seq, std::tie(view1, view2, view3), function);
  };

  const auto forEachLevel = [&](const auto& operation) {
    const auto function = [&](auto&& slice1, auto&& slice2, auto&& slice3) {
      const idx_t size = slice1.shape(0);
      for (idx_t idx = 0; idx < size; ++idx) {
        operation(slice1(idx), slice2(idx), slice3(idx));
      }
    };
    ArrayForEach<1>::apply(execution::seq, std::tie(view1, view2, view3), function);
  };

  const auto forEachAll = [&](const auto& operation) {
    ArrayForEach<0, 1>::apply(execution::seq, std::tie(view1, view2, view3), operation);
  };

  const auto forEachNested = [&](const auto& operation) {
      const auto function = [&](auto&& slice1, auto&& slice2, auto&& slice3) {
          ArrayForEach<0>::apply(execution::seq, std::tie(slice1, slice2, slice3), operation);
      };
      ArrayForEach<0>::apply(execution::seq, std::tie(view1, view2, view3), function);
  };

  const auto forEachConf = [&](const auto& operation) {
      const auto function = [&](auto&& slice1, auto&& slice2, auto&& slice3) {
          ArrayForEach<0>::apply(option::execution_policy(execution::seq), std::tie(slice1, slice2, slice3), operation);
      };
      ArrayForEach<0>::apply(option::execution_policy(execution::seq), std::tie(view1, view2, view3), function);
  };

  double baseline;
  baseline = timeLoop(rawPointer, num_iter, num_first, add, 0, "Addition; raw pointer               ");
  timeLoop(ijLoop, num_iter, num_first, add, baseline, "Addition; for loop (i, j)           ");
  timeLoop(jiLoop, num_iter, num_first, add, baseline, "Addition; for loop (j, i)           ");
  timeLoop(forEachCol, num_iter, num_first, add, baseline, "Addition; for each (columns)        ");
  timeLoop(forEachLevel, num_iter, num_first, add, baseline, "Addition; for each (levels)         ");
  timeLoop(forEachAll, num_iter, num_first, add, baseline, "Addition; for each (all elements)   ");
  timeLoop(forEachNested, num_iter, num_first, add, baseline, "Addition; for each (nested)         ");
  timeLoop(forEachConf, num_iter, num_first, add, baseline, "Addition; for each (nested, config) ");
  Log::info() << std::endl;

  num_first = 2;
  num_iter = 5;
  baseline = timeLoop(rawPointer, num_iter, num_first, trig, 0, "Trig    ; raw pointer               ");
  timeLoop(ijLoop, num_iter, num_first, trig, baseline, "Trig    ; for loop (i, j)           ");
  timeLoop(jiLoop, num_iter, num_first, trig, baseline, "Trig    ; for loop (j, i)           ");
  timeLoop(forEachCol, num_iter, num_first, trig, baseline, "Trig    ; for each (columns)        ");
  timeLoop(forEachLevel, num_iter, num_first, trig, baseline, "Trig    ; for each (levels)         ");
  timeLoop(forEachAll, num_iter, num_first, trig, baseline, "Trig    ; for each (all elements)   ");
  timeLoop(forEachNested, num_iter, num_first, trig, baseline, "Trig    ; for each (nested)         ");
  timeLoop(forEachConf, num_iter, num_first, trig, baseline, "Trig    ; for each (nested, config) ");
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
