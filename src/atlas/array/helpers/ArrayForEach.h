/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <tuple>
#include <type_traits>
#include <string_view>

#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/Range.h"
#include "atlas/array/helpers/ArraySlicer.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {

namespace execution {

// As in C++17 std::execution namespace. Note: unsequenced_policy is a C++20 addition.
class sequenced_policy {};
class unsequenced_policy {};
class parallel_unsequenced_policy {};
class parallel_policy {};

// execution policy objects as in C++ std::execution namespace. Note: unseq is a C++20 addition.
inline constexpr sequenced_policy            seq{ /*unspecified*/ };
inline constexpr parallel_policy             par{ /*unspecified*/ };
inline constexpr parallel_unsequenced_policy par_unseq{ /*unspecified*/ };
inline constexpr unsequenced_policy          unseq{ /*unspecified*/ };


// Type names for execution policy (Not in C++ standard)
template <typename execution_policy>
constexpr std::string_view policy_name() {
  return "unsupported";
};
template <>
constexpr std::string_view policy_name<sequenced_policy>() {
  return "sequenced_policy";
};
template <>
constexpr std::string_view policy_name<unsequenced_policy>() {
  return "unsequenced_policy";
};
template <>
constexpr std::string_view policy_name<parallel_unsequenced_policy>() {
  return "parallel_unsequenced_policy";
};

template <typename execution_policy>
constexpr std::string_view policy_name(execution_policy) {
  return policy_name<execution_policy>();
}

// Type check for execution policy (Not in C++ standard)
template <typename execution_policy>
using is_execution_policy = std::enable_if_t<
    std::is_same_v<execution_policy, sequenced_policy> ||
    std::is_same_v<execution_policy, parallel_policy> ||
    std::is_same_v<execution_policy, parallel_unsequenced_policy> ||
    std::is_same_v<execution_policy, unsequenced_policy>>*;

}  // namespace execution

namespace option {

// Convert execution_policy objects to a util::Config

template <typename T>
util::Config execution_policy() {
  return util::Config("execution_policy", execution::policy_name<T>());
}

template <typename T>
util::Config execution_policy(T) {
  return execution_policy<T>();
}

} // namespace option

namespace array {
namespace helpers {

namespace detail {

template <typename>
struct IsTupleImpl : std::false_type {};

template <typename... Args>
struct IsTupleImpl<std::tuple<Args...>> : std::true_type {};

template <typename Tuple>
using IsTuple = std::enable_if_t<IsTupleImpl<Tuple>::value>*;

template <typename... Ts, typename T>
constexpr auto tuplePushBack(const std::tuple<Ts...>& tuple, T value) {
  return std::tuple_cat(tuple, std::make_tuple(value));
}

template <typename Policy, typename Functor>
void forEach(idx_t idxMax, const Functor& functor) {

  if constexpr(std::is_same_v<Policy, execution::parallel_unsequenced_policy>) {
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

template<int Rank>
struct GetRankImpl {
  constexpr static int rank = Rank;
};

template<typename>
struct GetRank;

template<template <typename, int> typename View, typename Value, int Rank>
struct GetRank<View<Value, Rank>> : GetRankImpl<Rank> {};
template<template <typename, int> typename View, typename Value, int Rank>
struct GetRank<View<Value, Rank>&> : GetRankImpl<Rank> {};
template<template <typename, int> typename View, typename Value, int Rank>
struct GetRank<const View<Value, Rank>&> : GetRankImpl<Rank> {};
template<template <typename, int> typename View, typename Value, int Rank>
struct GetRank<View<Value, Rank>&&> : GetRankImpl<Rank> {};
template<template <typename, int> typename View, typename Value, int Rank>
struct GetRank<const View<Value, Rank>&&> : GetRankImpl<Rank> {};

template <size_t ViewIdx = 0, typename... SlicerArgs, typename ArrayViewTuple>
auto makeSlices(const std::tuple<SlicerArgs...>& slicerArgs,
                ArrayViewTuple&& arrayViews) {

  using ViewTupleType = std::decay_t<ArrayViewTuple>;

  if constexpr(ViewIdx < std::tuple_size_v<ViewTupleType>) {

    auto&& view = std::get<ViewIdx>(arrayViews);
    using View = std::tuple_element_t<ViewIdx, ViewTupleType>;

    constexpr auto Dim = sizeof...(SlicerArgs);
    constexpr auto Rank = GetRank<View>::rank;
    const auto paddedArgs =
      std::tuple_cat(slicerArgs, argPadding<Rank - Dim>());

    const auto slicer = [&view](const auto&... args) {
      return std::make_tuple(view.slice(args...));
    };

    return std::tuple_cat(std::apply(slicer, paddedArgs),
      makeSlices<ViewIdx + 1>(slicerArgs, std::forward<ArrayViewTuple>(arrayViews)));
  }
  else {
    return std::make_tuple();
  }
}

template <typename ExecutionPolicy, int Dim, int... ItrDims>
struct ArrayForEachImpl;

template <typename ExecutionPolicy, int Dim, int ItrDim, int... ItrDims>
struct ArrayForEachImpl<ExecutionPolicy, Dim, ItrDim, ItrDims...> {
  template <typename ArrayViewTuple, typename Mask, typename Function,
            typename... SlicerArgs, typename... MaskArgs>
  static void apply(ArrayViewTuple&& arrayViews,
                    const Mask& mask,
                    const Function& function,
                    const std::tuple<SlicerArgs...>& slicerArgs,
                    const std::tuple<MaskArgs...>& maskArgs) {

    using namespace detail;

    // Iterate over this dimension.
    if constexpr(Dim == ItrDim) {

      // Get size of iteration dimenion from first view argument.
      const auto idxMax = std::get<0>(arrayViews).shape(ItrDim);

      forEach<ExecutionPolicy>(idxMax, [&](idx_t idx) {

        // Decay from parallel_unsequenced to unsequenced policy
        if constexpr(std::is_same_v<ExecutionPolicy, execution::parallel_unsequenced_policy>) {
          ArrayForEachImpl<execution::unsequenced_policy, Dim + 1, ItrDims...>::apply(
              std::forward<ArrayViewTuple>(arrayViews), mask, function,
              tuplePushBack(slicerArgs, idx),
              tuplePushBack(maskArgs, idx));
        }
        else {
          // Retain current execution policy.
          ArrayForEachImpl<ExecutionPolicy, Dim + 1, ItrDims...>::apply(
              std::forward<ArrayViewTuple>(arrayViews), mask, function,
              tuplePushBack(slicerArgs, idx),
              tuplePushBack(maskArgs, idx));
        }
      });
    }
    // Add a RangeAll to arguments.
    else {
      ArrayForEachImpl<ExecutionPolicy, Dim + 1, ItrDim, ItrDims...>::apply(
          std::forward<ArrayViewTuple>(arrayViews), mask, function,
          tuplePushBack(slicerArgs, Range::all()),
          maskArgs);
    }
  }
};

template <typename ExecutionPolicy, int Dim>
struct ArrayForEachImpl<ExecutionPolicy, Dim> {
  template <typename ArrayViewTuple, typename Mask, typename Function,
            typename... SlicerArgs, typename... MaskArgs>
  static void apply(ArrayViewTuple&& arrayViews,
                    const Mask& mask,
                    const Function& function,
                    const std::tuple<SlicerArgs...>& slicerArgs,
                    const std::tuple<MaskArgs...>& maskArgs) {

    // Skip iteration if mask evaluates to true.
    if (std::apply(mask, maskArgs)) {
      return;
    }

    auto slices = makeSlices(slicerArgs, std::forward<ArrayViewTuple>(arrayViews));
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
  /// details Visits all elements indexed by ItrDims and creates a slice from
  ///         each ArrayView in arrayViews. Slices are sent to function
  ///         which is executed with the signature f(slice1, slice2,...).
  ///         Iterations are skipped when mask evaluates to "true"
  ///         and is executed with signature g(idx_i, idx_j,...), where the idxs
  ///         are indices of ItrDims.
  ///         When a config is supplied containing "execution_policy" =
  ///         "parallel_unsequenced_policy" (default) the first loop is executed
  ///         using OpenMP. The remaining loops are executed in serial. When
  ///         "execution_policy" = "sequenced_policy", all loops are executed in
  ///         sequential (row-major) order.
  ///         Note: The lowest ArrayView.rank() must be greater than or equal
  ///         to the highest dim in ItrDims. TODO: static checking for this.
  template <typename ArrayViewTuple, typename Mask, typename Function,
            detail::IsTuple<ArrayViewTuple> = nullptr>
  static void apply(const eckit::Parametrisation& conf,
                    ArrayViewTuple&& arrayViews,
                    const Mask& mask, const Function& function) {

    auto execute = [&](auto execution_policy) {
      apply(execution_policy, std::forward<ArrayViewTuple>(arrayViews), mask, function);
    };

    using namespace execution;
    std::string execution_policy;
    if( conf.get("execution_policy",execution_policy) ) {
      if (execution_policy == policy_name(par_unseq)) {
        execute(par_unseq);
      } else if (execution_policy == policy_name(par)) {
        execute(par);
      } else if (execution_policy == policy_name(unseq)) {
        execute(unseq);
      } else if (execution_policy == policy_name(seq)) {
        execute(seq);
      }
    }
    else {
      execute(par_unseq);
    }
  }

  /// brief   Apply "For-Each" method.
  ///
  /// details As above, but Execution policy is determined at compile-time.
  template <typename ExecutionPolicy, typename ArrayViewTuple, typename Mask, typename Function,
            execution::is_execution_policy<ExecutionPolicy> = nullptr,
            detail::IsTuple<ArrayViewTuple> = nullptr>
  static void apply(ExecutionPolicy, ArrayViewTuple&& arrayViews, const Mask& mask, const Function& function) {

    detail::ArrayForEachImpl<ExecutionPolicy, 0, ItrDims...>::apply(
        std::forward<ArrayViewTuple>(arrayViews), mask, function, std::make_tuple(), std::make_tuple());
  }

  /// brief   Apply "For-Each" method
  ///
  /// detials Apply ForEach with default execution policy.
  template <typename ArrayViewTuple, typename Mask, typename Function,
            detail::IsTuple<ArrayViewTuple> = nullptr>
  static void apply(ArrayViewTuple&& arrayViews, const Mask& mask, const Function& function) {
      apply(std::forward<ArrayViewTuple>(arrayViews), mask, function);
  }

  /// brief   Apply "For-Each" method
  ///
  /// detials Apply ForEach with run-time determined execution policy and no mask.
  template <typename ArrayViewTuple, typename Function,
            detail::IsTuple<ArrayViewTuple> = nullptr>
  static void apply(const eckit::Parametrisation& conf, ArrayViewTuple&& arrayViews, const Function& function) {
    constexpr auto no_mask = [](auto args...) { return 0; };
    apply(conf, std::forward<ArrayViewTuple>(arrayViews), no_mask, function);
  }

  /// brief   Apply "For-Each" method
  ///
  /// detials Apply ForEach with compile-time determined execution policy and no mask.
  template <typename ExecutionPolicy, typename ArrayViewTuple,
            typename Function, detail::IsTuple<ArrayViewTuple> = nullptr,
            execution::is_execution_policy<ExecutionPolicy> = nullptr>
  static void apply(ExecutionPolicy executionPolicy, ArrayViewTuple&& arrayViews, const Function& function) {
    constexpr auto no_mask = [](auto args...) { return 0; };
    apply(executionPolicy, std::forward<ArrayViewTuple>(arrayViews), no_mask, function);
  }

  /// brief   Apply "For-Each" method
  ///
  /// detials Apply ForEach with default execution policy and no mask.
  template <typename ArrayViewTuple, typename Function,
            detail::IsTuple<ArrayViewTuple> = nullptr>
  static void apply(ArrayViewTuple arrayViews, const Function& function) {
    apply(execution::par_unseq, std::forward<ArrayViewTuple>(arrayViews), function);
  }

};

}  // namespace helpers
}  // namespace array
}  // namespace atlas
