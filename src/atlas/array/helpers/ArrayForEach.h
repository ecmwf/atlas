/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>
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

namespace detail {

enum class ExecutionPolicy {
  parallel,
  serial
};

template <typename... Ts, typename T>
constexpr auto tuplePushBack(const std::tuple<Ts...>& tuple, T value) {
  return std::tuple_cat(tuple, std::make_tuple(value));
}

template <ExecutionPolicy Policy, typename Functor>
void forEach(idx_t idxMax, const Functor& functor) {

  if constexpr(Policy == ExecutionPolicy::parallel) {
      atlas_omp_parallel_for(auto idx = idx_t{}; idx < idxMax; ++idx) {
        functor(idx);
      }
    }
  else {
    for (auto idx = idx_t{}; idx < idxMax; ++idx) {
      functor(idx);
    }
  }
}

template <size_t ViewIdx = 0, typename... ArrayViews, typename... SlicerArgs>
auto makeSlices(std::tuple<ArrayViews...>& arrayViews,
                const std::tuple<SlicerArgs...>& slicerArgs) {

  // Loop over all views in tuple and make slices.
  if constexpr(ViewIdx < std::tuple_size_v<std::tuple<ArrayViews...>>) {

      auto& view = std::get<ViewIdx>(arrayViews);
      constexpr auto rank = std::get<ViewIdx>(arrayViews).rank();

      // "Lambdafy" slicer apply method to work with std::apply.
      const auto slicer = [&view](const auto&... args) {
        return view.slice(args...);
      };

      // Fill out the remaining slicerArgs with Range::all().
      constexpr auto dim = std::tuple_size_v<std::tuple<SlicerArgs...>>;
      auto argPadding = std::array<decltype(Range::all()), rank - dim>();
      argPadding.fill(Range::all());
      const auto paddedArgs = std::tuple_cat(slicerArgs, argPadding);

      // Recurse until all views are sliced.
      return std::tuple_cat(std::make_tuple(std::apply(slicer, paddedArgs)),
                            makeSlices<ViewIdx + 1>(arrayViews, slicerArgs));
    }
  else {

    // No more ArrayViews. Finish recursion.
    return std::make_tuple();
  }
}

template <ExecutionPolicy Policy, int Dim, int... ItrDims>
struct ArrayForEachImpl;

template <ExecutionPolicy Policy, int Dim, int ItrDim, int... ItrDims>
struct ArrayForEachImpl<Policy, Dim, ItrDim, ItrDims...> {
  template <typename... ArrayViews, typename LoopFunctor, typename GhostFunctor,
            typename... SlicerArgs, typename... GhostArgs>
  static void apply(std::tuple<ArrayViews...>& arrayViews,
                    const LoopFunctor& loopFunctor,
                    const GhostFunctor& ghostFunctor,
                    const std::tuple<SlicerArgs...>& slicerArgs,
                    const std::tuple<GhostArgs...>& ghostArgs) {

    using namespace detail;

    // Iterate over this dimension.
    if constexpr(Dim == ItrDim) {

        // Get size of iteration dimenion from first view argument.
        const auto idxMax = std::get<0>(arrayViews).shape(ItrDim);

        forEach<Policy>(idxMax, [&](idx_t idx) {

          // Always set Policy to serial after one loop.
          ArrayForEachImpl<ExecutionPolicy::serial, Dim + 1,
                           ItrDims...>::apply(arrayViews, loopFunctor,
                                              ghostFunctor,
                                              tuplePushBack(slicerArgs, idx),
                                              tuplePushBack(ghostArgs, idx));
        });
      }

    // Add a RangeAll to arguments.
    else {
      ArrayForEachImpl<Policy, Dim + 1, ItrDim, ItrDims...>::apply(
          arrayViews, loopFunctor, ghostFunctor,
          tuplePushBack(slicerArgs, Range::all()), ghostArgs);
    }
  }
};

template <ExecutionPolicy Policy, int Dim>
struct ArrayForEachImpl<Policy, Dim> {
  template <typename... ArrayViews, typename LoopFunctor, typename GhostFunctor,
            typename... SlicerArgs, typename... GhostArgs>
  static void apply(std::tuple<ArrayViews...>& arrayViews,
                    const LoopFunctor& loopFunctor,
                    const GhostFunctor& ghostFunctor,
                    const std::tuple<SlicerArgs...>& slicerArgs,
                    const std::tuple<GhostArgs...>& ghostArgs) {

    // Skip iteration if ghostFunctor evaluates to true.
    if (std::apply(ghostFunctor, ghostArgs)) {
      return;
    }

    auto slices = makeSlices(arrayViews, slicerArgs);
    std::apply(loopFunctor, slices);
  }
};
}  // namespace detail


/// brief  Array "For-Each" helper struct.
///
/// detail Iterates over dimensions given in  ItrDims. Slices over full range
/// of other dimensions.
/// Note: ItrDims must be given in ascending numerical order. TODO: Static
/// checking for this.
template <int... ItrDims>
struct ArrayForEach {
  /// brief   Apply "For-Each" method.
  ///
  /// detials Vists all elements indexed by ItrDims and creates a slice from
  ///         each ArrayView in arrayViews. Slices are sent to loopFunctor
  ///         which is executed with the signature f(slice1, slice2,...).
  ///         Iterations are skipped when ghostFunctor evaluates to "true"
  ///         and is executed with signature g(idx_i, idx_j,...), where the idxs
  ///         are indices of ItrDims.
  ///         When the config contains "execution policy" = "parallel" (default)
  ///         the first loop is executed using OpenMP. The remaining loops are
  ///         executed in serial. When "execution policy" = "serial", all loops
  ///         are executed in serial.
  ///         Note: The lowest ArrayView.rank() must be greater than or equal
  ///         to the highest dim in ItrDims. TODO: static checking for this.
  template <typename... ArrayViews, typename LoopFunctor, typename GhostFunctor>
  static void apply(const std::tuple<ArrayViews...>& arrayViews,
                    const LoopFunctor& loopFunctor,
                    const GhostFunctor& ghostFunctor,
                    const util::Config& conf = util::Config()) {

    using namespace detail;

    // Make a copy of views to simplify constness and forwarding.
    auto arrayViewsCopy = arrayViews;

    const auto executionPolicy = conf.getString("execution policy", "parallel");

    if (executionPolicy == "parallel") {
      ArrayForEachImpl<ExecutionPolicy::parallel, 0, ItrDims...>::apply(
          arrayViewsCopy, loopFunctor, ghostFunctor, std::make_tuple(),
          std::make_tuple());
    } else if (executionPolicy == "serial") {
      ArrayForEachImpl<ExecutionPolicy::serial, 0, ItrDims...>::apply(
          arrayViewsCopy, loopFunctor, ghostFunctor, std::make_tuple(),
          std::make_tuple());
    } else {
      throw eckit::BadParameter("Unknown execution policy: " + executionPolicy,
                                Here());
    }
  }

  /// brief   Apply "For-Each" method.
  ///
  /// details Apply "For-Each" without a ghostFunctor.
  template <typename... ArrayViews, typename LoopFunctor>
  static void apply(const std::tuple<ArrayViews...>& arrayViews,
                    const LoopFunctor& loopFunctor,
                    const util::Config& conf = util::Config()) {
    apply(arrayViews, loopFunctor, [](auto args...) { return 0; }, conf);
  }
};

}  // namespace helpers
}  // namespace array
}  // namespace atlas
