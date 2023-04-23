/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <tuple>

#include "atlas/array/ArrayView.h"
#include "atlas/array/Range.h"
#include "atlas/array/helpers/ArraySlicer.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/util/Config.h"

#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace array {
namespace helpers {

enum class execution_policy {
  parallel,
  serial
};

namespace detail {

template <typename... Ts, typename T>
constexpr auto tupleAppend(const std::tuple<Ts...>& tuple, const T value) {
  return std::tuple_cat(tuple, std::make_tuple(value));
}

template <execution_policy Policy, typename Functor>
void forEach(idx_t idxMax, const Functor& functor) {

    if constexpr (Policy == execution_policy::parallel) {
        atlas_omp_parallel_for (auto idx = idx_t{}; idx < idxMax; ++idx) {
            functor(idx);
        }
    } else {
        for (auto idx = idx_t{}; idx < idxMax; ++idx) {
            functor(idx);
        }
    }
}

template <int N>
constexpr auto rangeAlls() {
  return std::tuple_cat(std::make_tuple(Range::all()), rangeAlls<N - 1>());
}

template <>
constexpr auto rangeAlls<0>() {
  return std::make_tuple();
};

template <int Dim, typename Value, int Rank, typename... SlicerArgs>
auto makeSlice(ArrayView<Value, Rank>& view,
                const std::tuple<SlicerArgs...>& slicerArgs) {

  using View = typename std::remove_reference<decltype(view)>::type;

  // "Lambdafy" slicer apply method to work with std::apply.
  const auto slicer = [&view](const auto&... args) {
    return ArraySlicer<View>(view).apply(args...);
  };

  const auto paddedArgs =
      std::tuple_cat(slicerArgs, rangeAlls<Rank - Dim>());

  return std::apply(slicer, paddedArgs);
}

}  // namespace detail

template <execution_policy Policy, int Dim, int... ItrDims>
struct ArrayForEachImpl;

template <execution_policy Policy, int Dim, int ItrDim, int... ItrDims>
struct ArrayForEachImpl<Policy, Dim, ItrDim, ItrDims...> {
  template <typename Value1, typename Value2, int Rank1, int Rank2,
            typename LoopFunctor, typename GhostFunctor,
            typename... SlicerArgs, typename... GhostArgs>
  static void apply(ArrayView<Value1, Rank1>& view1,
                    ArrayView<Value2, Rank2>& view2,
                    const LoopFunctor& loopFunctor,
                    const GhostFunctor& ghostFunctor,
                    const std::tuple<SlicerArgs...>& slicerArgs,
                    const std::tuple<GhostArgs...>& ghostArgs) {

    // Iterate over this dimension.
    if constexpr(Dim == ItrDim) {
        detail::forEach<Policy>(view1.shape(ItrDim), [&](idx_t idx){
            ArrayForEachImpl<execution_policy::serial, Dim + 1, ItrDims...>::
                apply(view1, view2, loopFunctor, ghostFunctor,
                detail::tupleAppend(slicerArgs, idx), detail::tupleAppend(ghostArgs, idx));
        });
      }

    // Add aa RangeAll to arguments and skip.
    else {
      ArrayForEachImpl<Policy, Dim + 1, ItrDim, ItrDims...>::apply(
          view1, view2, loopFunctor, ghostFunctor,
          detail::tupleAppend(slicerArgs, Range::all()), ghostArgs);
    }
  }
};

template <execution_policy Policy, int Dim>
struct ArrayForEachImpl<Policy, Dim> {
  template <typename Value1, typename Value2, int Rank1, int Rank2,
            typename LoopFunctor, typename GhostFunctor,
            typename... SlicerArgs, typename... GhostArgs>
  static void apply(ArrayView<Value1, Rank1>& view1,
                      ArrayView<Value2, Rank2>& view2,
                      const LoopFunctor& loopFunctor,
                      const GhostFunctor& ghostFunctor,
                      const std::tuple<SlicerArgs...>& slicerArgs,
                      const std::tuple<GhostArgs...>& ghostArgs) {

    // Don't apply loopFunctor if ghostFunctor returns true.
    if (std::apply(ghostFunctor, ghostArgs)) {
      return;
    }

    auto slice1 = detail::makeSlice<Dim>(view1, slicerArgs);
    auto slice2 = detail::makeSlice<Dim>(view2, slicerArgs);

    loopFunctor(slice1, slice2);
  }
};

template <int... ItrDims>
struct ArrayForEach {
  template <typename Value1, typename Value2, int Rank1, int Rank2,
            typename LoopFunctor, typename GhostFunctor>
  static void apply(ArrayView<Value1, Rank1> view1,
                    ArrayView<Value2, Rank2> view2,
                    const LoopFunctor& loopFunctor,
                    const GhostFunctor& ghostFunctor,
                    const util::Config& conf = util::Config()) {

    const auto executionPolicy = conf.getString("execution policy", "parallel");

    if (executionPolicy == "parallel") {
        ArrayForEachImpl<execution_policy::parallel, 0, ItrDims...>::apply(
            view1, view2, loopFunctor, ghostFunctor, std::make_tuple(), std::make_tuple());
    } else if (executionPolicy == "serial") {
        ArrayForEachImpl<execution_policy::serial, 0, ItrDims...>::apply(
            view1, view2, loopFunctor, ghostFunctor, std::make_tuple(), std::make_tuple());
    } else {
        throw eckit::BadParameter("Unknown execution policy: " + executionPolicy, Here());
    }
  }

  template <typename Value1, typename Value2, int Rank1, int Rank2, typename LoopFunctor>
  static void apply(ArrayView<Value1, Rank1> view1,
                    ArrayView<Value2, Rank2> view2,
                    const LoopFunctor& loopFunctor,
                    const util::Config& conf = util::Config()) {
      apply(view1, view2, loopFunctor, [](const auto args){return 0;});
  }
};

}  // namespace helpers
}  // namespace array
}  // namespace atlas
