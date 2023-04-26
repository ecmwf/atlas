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

auto sequencedConf() {
  return util::Config("execution_policy", "sequenced");
}

auto unsequencedConf() {
  return util::Config("execution_policy", "unsequenced");
}

auto parallelUnsequencedConf() {
  return util::Config("execution_policy", "parallel_unsequenced");
}

enum class ExecutionPolicy {
  parallel_unsequenced_policy,
  unsequenced_policy,
  sequenced_policy
};

template <typename... Ts, typename T>
constexpr auto tuplePushBack(const std::tuple<Ts...>& tuple, T value) {
  return std::tuple_cat(tuple, std::make_tuple(value));
}

template <ExecutionPolicy Policy, typename Functor>
void forEach(idx_t idxMax, const Functor& functor) {

  if constexpr(Policy == ExecutionPolicy::parallel_unsequenced_policy) {
      atlas_omp_parallel_for(auto idx = idx_t{}; idx < idxMax; ++idx) {
        functor(idx);
      }
    }
  else {
    // Simple for-loop for sequenced or unsequenced execution policies.
    for (auto idx = idx_t{}; idx < idxMax; ++idx) {
      functor(idx);
    }
  }
}

template <int NPad>
constexpr auto argPadding() {
  if constexpr(NPad > 0) {
      return std::tuple_cat(std::make_tuple(Range::all()),
                            argPadding<NPad - 1>());
    }
  else {
    return std::make_tuple();
  }
}

template <typename... SlicerArgs, template <typename, int> typename View,
          typename Value, int Rank, typename... ArrayViews>
auto makeSlices(const std::tuple<SlicerArgs...>& slicerArgs,
                View<Value, Rank>& arrayView, ArrayViews&... arrayViews) {

  // "Lambdafy" slicer apply method to work with std::apply.
  const auto slicer = [&arrayView](const auto&... args) {
    return arrayView.slice(args...);
  };

  // Fill out the remaining slicerArgs with Range::all().
  constexpr auto Dim = sizeof...(SlicerArgs);
  const auto paddedArgs = std::tuple_cat(slicerArgs, argPadding<Rank - Dim>());

  if constexpr(sizeof...(ArrayViews) > 0) {

      // Recurse until all views are sliced.
      return std::tuple_cat(std::make_tuple(std::apply(slicer, paddedArgs)),
                            makeSlices(slicerArgs, arrayViews...));
    }
  else {
    return std::make_tuple(std::apply(slicer, paddedArgs));
  }
}

template <ExecutionPolicy Policy, int Dim, int... ItrDims>
struct ArrayForEachImpl;

template <ExecutionPolicy Policy, int Dim, int ItrDim, int... ItrDims>
struct ArrayForEachImpl<Policy, Dim, ItrDim, ItrDims...> {
  template <typename... ArrayViews, typename Function, typename Mask,
            typename... SlicerArgs, typename... MaskArgs>
  static void apply(std::tuple<ArrayViews...>& arrayViews,
                    const Function& function, const Mask& mask,
                    const std::tuple<SlicerArgs...>& slicerArgs,
                    const std::tuple<MaskArgs...>& maskArgs) {

    using namespace detail;

    // Iterate over this dimension.
    if constexpr(Dim == ItrDim) {

        // Get size of iteration dimenion from first view argument.
        const auto idxMax = std::get<0>(arrayViews).shape(ItrDim);

        forEach<Policy>(idxMax, [&](idx_t idx) {

          // Decay from parallel_unsequenced to unsequenced policy
          if constexpr(Policy == ExecutionPolicy::parallel_unsequenced_policy) {
              ArrayForEachImpl<ExecutionPolicy::unsequenced_policy, Dim + 1,
                               ItrDims...>::apply(arrayViews, function, mask,
                                                  tuplePushBack(slicerArgs,
                                                                idx),
                                                  tuplePushBack(maskArgs, idx));
            }
          else {
            // Retain current execution policy.
            ArrayForEachImpl<Policy, Dim + 1, ItrDims...>::apply(
                arrayViews, function, mask, tuplePushBack(slicerArgs, idx),
                tuplePushBack(maskArgs, idx));
          }
        });
      }

    // Add a RangeAll to arguments.
    else {
      ArrayForEachImpl<Policy, Dim + 1, ItrDim, ItrDims...>::apply(
          arrayViews, function, mask, tuplePushBack(slicerArgs, Range::all()),
          maskArgs);
    }
  }
};

template <ExecutionPolicy Policy, int Dim>
struct ArrayForEachImpl<Policy, Dim> {
  template <typename... ArrayViews, typename Function, typename Mask,
            typename... SlicerArgs, typename... MaskArgs>
  static void apply(std::tuple<ArrayViews...>& arrayViews,
                    const Function& function, const Mask& mask,
                    const std::tuple<SlicerArgs...>& slicerArgs,
                    const std::tuple<MaskArgs...>& maskArgs) {

    // Skip iteration if mask evaluates to true.
    if (std::apply(mask, maskArgs)) {
      return;
    }

    const auto slicerWrapper = [&slicerArgs](auto&... args) {
      return makeSlices(slicerArgs, args...);
    };

    auto slices = std::apply(slicerWrapper, arrayViews);
    std::apply(function, slices);
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
  ///         each ArrayView in arrayViews. Slices are sent to function
  ///         which is executed with the signature f(slice1, slice2,...).
  ///         Iterations are skipped when mask evaluates to "true"
  ///         and is executed with signature g(idx_i, idx_j,...), where the idxs
  ///         are indices of ItrDims.
  ///         When the config contains "execution_policy" =
  ///         "parallel_unsequenced" (default) the first loop is executed using
  ///         OpenMP. The remaining loops are executed in serial. When
  ///         "execution_policy" = "sequenced", all loops are executed in
  ///         sequential (row-major) order.
  ///         Note: The lowest ArrayView.rank() must be greater than or equal
  ///         to the highest dim in ItrDims. TODO: static checking for this.
  template <typename... ArrayViews, typename Function, typename Mask>
  static void apply(const std::tuple<ArrayViews...>& arrayViews,
                    const Function& function, const Mask& mask,
                    const util::Config& conf =
                        detail::parallelUnsequencedConf()) {

    using namespace detail;

    // Make a copy of views to simplify constness and forwarding.
    auto arrayViewsCopy = arrayViews;

    const auto executionPolicy = conf.getString("execution_policy");

    if (executionPolicy == "parallel_unsequenced") {
      ArrayForEachImpl<ExecutionPolicy::parallel_unsequenced_policy, 0,
                       ItrDims...>::apply(arrayViewsCopy, function, mask,
                                          std::make_tuple(), std::make_tuple());
    } else if (executionPolicy == "unsequenced") {
      ArrayForEachImpl<ExecutionPolicy::unsequenced_policy, 0,
                       ItrDims...>::apply(arrayViewsCopy, function, mask,
                                          std::make_tuple(), std::make_tuple());
    } else if (executionPolicy == "sequenced") {
      ArrayForEachImpl<ExecutionPolicy::sequenced_policy, 0, ItrDims...>::apply(
          arrayViewsCopy, function, mask, std::make_tuple(), std::make_tuple());
    } else {
      throw eckit::BadParameter("Unknown execution policy: " + executionPolicy,
                                Here());
    }
  }

  /// brief   Apply "For-Each" method.
  ///
  /// details Apply "For-Each" without a mask.
  template <typename... ArrayViews, typename Function>
  static void apply(const std::tuple<ArrayViews...>& arrayViews,
                    const Function& function,
                    const util::Config& conf =
                        detail::parallelUnsequencedConf()) {
    apply(arrayViews, function, [](auto args...) { return 0; }, conf);
  }
};

}  // namespace helpers
}  // namespace array
}  // namespace atlas
