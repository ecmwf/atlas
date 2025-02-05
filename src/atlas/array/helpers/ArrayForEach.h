/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

#include "atlas/array/ArrayView.h"
#include "atlas/array/Range.h"
#include "atlas/array/helpers/ArraySlicer.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {

namespace execution {

// As in C++17 std::execution namespace. Note: unsequenced_policy is a C++20
// addition.
class sequenced_policy {};
class unsequenced_policy {};
class parallel_unsequenced_policy {};
class parallel_policy {};

// execution policy objects as in C++ std::execution namespace. Note: unseq is a
// C++20 addition.
inline constexpr sequenced_policy seq{/*unspecified*/};
inline constexpr parallel_policy par{/*unspecified*/};
inline constexpr parallel_unsequenced_policy par_unseq{/*unspecified*/};
inline constexpr unsequenced_policy unseq{/*unspecified*/};

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
template <>
constexpr std::string_view policy_name<parallel_policy>() {
  return "parallel_policy";
};

template <typename execution_policy>
constexpr std::string_view policy_name(execution_policy) {
  return policy_name<std::decay_t<execution_policy>>();
}

// Type check for execution policy (Not in C++ standard)
template <typename execution_policy>
constexpr auto is_execution_policy() {
  return std::is_same_v<execution_policy, sequenced_policy> ||
         std::is_same_v<execution_policy, parallel_policy> ||
         std::is_same_v<execution_policy, parallel_unsequenced_policy> ||
         std::is_same_v<execution_policy, unsequenced_policy>;
}

template <typename ExecutionPolicy>
constexpr auto demote_policy() {
  if constexpr (std::is_same_v<ExecutionPolicy, parallel_unsequenced_policy>) {
    return unseq;
  } else if constexpr (std::is_same_v<ExecutionPolicy, parallel_policy>) {
    return seq;
  } else {
    return ExecutionPolicy{};
  }
  ATLAS_UNREACHABLE();
}

template <typename execution_policy>
constexpr auto is_omp_policy() {
  return std::is_same_v<execution_policy, parallel_policy> ||
         std::is_same_v<execution_policy, parallel_unsequenced_policy>;
}

template <typename ExecutionPolicy>
using demote_policy_t = decltype(demote_policy<ExecutionPolicy>());

}  // namespace execution

namespace option {

// Convert execution_policy objects to a util::Config

template <typename T>
util::Config execution_policy() {
  return util::Config("execution_policy",
                      execution::policy_name<std::decay_t<T>>());
}

template <typename T>
util::Config execution_policy(T) {
  return execution_policy<std::decay_t<T>>();
}

}  // namespace option

namespace array {
namespace helpers {

namespace detail {

struct NoMask {
  template <typename... Args>
  constexpr bool operator()(Args...) const {
    return 0;
  }
};

inline constexpr NoMask no_mask;

template <typename... Ts, typename T>
constexpr auto tuplePushBack(const std::tuple<Ts...>& tuple, T value) {
  return std::tuple_cat(tuple, std::make_tuple(value));
}

template <typename ExecutionPolicy, typename Functor>
void forEach(idx_t idxMax, const Functor& functor) {
  if constexpr (execution::is_omp_policy<ExecutionPolicy>()) {
    atlas_omp_parallel_for(auto idx = idx_t{}; idx < idxMax; ++idx) {
      functor(idx);
    }
  } else {
    // Simple for-loop for sequenced or unsequenced execution policies.
    for (auto idx = idx_t{}; idx < idxMax; ++idx) {
      functor(idx);
    }
  }
}

template <int NPad>
constexpr auto argPadding() {
  if constexpr (NPad > 0) {
    return std::tuple_cat(std::make_tuple(Range::all()),
                          argPadding<NPad - 1>());
  } else {
    return std::make_tuple();
  }
  ATLAS_UNREACHABLE();
}

template <size_t ViewIdx = 0, typename... SlicerArgs, typename ArrayViewTuple>
auto makeSlices(const std::tuple<SlicerArgs...>& slicerArgs,
                ArrayViewTuple&& arrayViews) {
  constexpr auto nb_views = std::tuple_size_v<ArrayViewTuple>;

  auto&& arrayView = std::get<ViewIdx>(arrayViews);
  using ArrayView = std::decay_t<decltype(arrayView)>;

  constexpr auto Dim = sizeof...(SlicerArgs);
  constexpr auto Rank = ArrayView::rank();

  static_assert(
      Dim <= Rank,
      "Error: number of slicer arguments exceeds the rank of ArrayView.");
  const auto paddedArgs = std::tuple_cat(slicerArgs, argPadding<Rank - Dim>());

  const auto slicer = [&arrayView](const auto&... args) {
    return std::make_tuple(arrayView.slice(args...));
  };

  if constexpr (ViewIdx == nb_views - 1) {
    return std::apply(slicer, paddedArgs);
  } else {
    // recurse
    return std::tuple_cat(
        std::apply(slicer, paddedArgs),
        makeSlices<ViewIdx + 1>(slicerArgs,
                                std::forward<ArrayViewTuple>(arrayViews)));
  }
  ATLAS_UNREACHABLE();
}

template <typename ExecutionPolicy, int Dim, int... ItrDims>
struct ArrayForEachImpl;

template <typename ExecutionPolicy, int Dim, int ItrDim, int... ItrDims>
struct ArrayForEachImpl<ExecutionPolicy, Dim, ItrDim, ItrDims...> {
  template <typename ArrayViewTuple, typename Mask, typename Function,
            typename... SlicerArgs, typename... MaskArgs>
  static void apply(ArrayViewTuple&& arrayViews, const Mask& mask,
                    const Function& function,
                    const std::tuple<SlicerArgs...>& slicerArgs,
                    const std::tuple<MaskArgs...>& maskArgs) {
    // Iterate over this dimension.
    if constexpr (Dim == ItrDim) {
      // Get size of iteration dimenion from first view argument.
      const auto idxMax = std::get<0>(arrayViews).shape(ItrDim);

      forEach<ExecutionPolicy>(idxMax, [&](idx_t idx) {
        // Demote parallel execution policy to a non-parallel one in further
        // recursion
        ArrayForEachImpl<
            execution::demote_policy_t<ExecutionPolicy>, Dim + 1,
            ItrDims...>::apply(std::forward<ArrayViewTuple>(arrayViews), mask,
                               function, tuplePushBack(slicerArgs, idx),
                               tuplePushBack(maskArgs, idx));
      });
    }
    // Add a RangeAll to arguments.
    else {
      ArrayForEachImpl<ExecutionPolicy, Dim + 1, ItrDim, ItrDims...>::apply(
          std::forward<ArrayViewTuple>(arrayViews), mask, function,
          tuplePushBack(slicerArgs, Range::all()), maskArgs);
    }
  }
};

template <typename...>
struct is_applicable : std::false_type {};

template <typename Function, typename... Args>
struct is_applicable<Function, std::tuple<Args...>>
    : std::is_invocable<Function, Args...> {};

template <typename Function, typename Tuple>
inline constexpr bool is_applicable_v = is_applicable<Function, Tuple>::value;

template <typename ExecutionPolicy, int Dim>
struct ArrayForEachImpl<ExecutionPolicy, Dim> {
  template <typename ArrayViewTuple, typename Mask, typename Function,
            typename... SlicerArgs, typename... MaskArgs>
  static void apply(ArrayViewTuple&& arrayViews, const Mask& mask,
                    const Function& function,
                    const std::tuple<SlicerArgs...>& slicerArgs,
                    const std::tuple<MaskArgs...>& maskArgs) {
    constexpr auto maskPresent = !std::is_same_v<Mask, NoMask>;

    if constexpr (maskPresent) {
      constexpr auto invocableMask =
          std::is_invocable_r_v<int, Mask, MaskArgs...>;
      static_assert(
          invocableMask,
          "Cannot invoke mask function with given arguments.\n"
          "Make sure you arguments are N integers (or auto...) "
          "where N == sizeof...(ItrDims). Function must return an int.");

      if (std::apply(mask, maskArgs)) {
        return;
      }
    }

    auto slices =
        makeSlices(slicerArgs, std::forward<ArrayViewTuple>(arrayViews));

    constexpr auto applicable = is_applicable_v<Function, decltype(slices)>;
    static_assert(
        applicable,
        "Cannot invoke function with given arguments. "
        "Make sure you the arguments are rvalue references (Slice&&) or const "
        "references (const Slice&) or regular value (Slice)");
    std::apply(function, std::move(slices));
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
  ///         "sequenced_policy" (default). All loops are then executed in
  ///         sequential (row-major) order. With "execution_policy" =
  ///         "parallel_unsequenced" the first loop is executed using OpenMP.
  ///         The remaining loops are executed in serial. Note: The lowest
  ///         ArrayView.rank() must be greater than or equal to the highest dim
  ///         in ItrDims. TODO: static checking for this.
  template <typename... ArrayView, typename Mask, typename Function>
  static void apply(const eckit::Parametrisation& conf,
                    std::tuple<ArrayView...>&& arrayViews, const Mask& mask,
                    const Function& function) {
    auto execute = [&](auto execution_policy) {
      apply(execution_policy, std::move(arrayViews), mask, function);
    };

    using namespace execution;
    std::string execution_policy;
    if (conf.get("execution_policy", execution_policy)) {
      if (execution_policy == policy_name(par_unseq)) {
        execute(par_unseq);
      } else if (execution_policy == policy_name(par)) {
        execute(par);
      } else if (execution_policy == policy_name(unseq)) {
        execute(unseq);
      } else if (execution_policy == policy_name(seq)) {
        execute(seq);
      } else {
        throw_Exception("Unrecognized execution policy " + execution_policy,
                        Here());
      }
    } else {
      execute(seq);
    }
  }

  /// brief   Apply "For-Each" method.
  ///
  /// details As above, but Execution policy is determined at compile-time.
  template <typename ExecutionPolicy, typename... ArrayView, typename Mask,
            typename Function,
            typename = std::enable_if_t<
                execution::is_execution_policy<ExecutionPolicy>()>>
  static void apply(ExecutionPolicy, std::tuple<ArrayView...>&& arrayViews,
                    const Mask& mask, const Function& function) {
    detail::ArrayForEachImpl<ExecutionPolicy, 0, ItrDims...>::apply(
        std::move(arrayViews), mask, function, std::make_tuple(),
        std::make_tuple());
  }

  /// brief   Apply "For-Each" method
  ///
  /// details Apply ForEach with default execution policy.
  template <typename... ArrayView, typename Mask, typename Function>
  static void apply(std::tuple<ArrayView...>&& arrayViews, const Mask& mask,
                    const Function& function) {
    apply(std::move(arrayViews), mask, function);
  }

  /// brief   Apply "For-Each" method
  ///
  /// details Apply ForEach with run-time determined execution policy and no
  /// mask.
  template <typename... ArrayView, typename Function>
  static void apply(const eckit::Parametrisation& conf,
                    std::tuple<ArrayView...>&& arrayViews,
                    const Function& function) {
    apply(conf, std::move(arrayViews), detail::no_mask, function);
  }

  /// brief   Apply "For-Each" method
  ///
  /// details Apply ForEach with compile-time determined execution policy and no
  /// mask.
  template <typename ExecutionPolicy, typename... ArrayView, typename Function,
            typename = std::enable_if_t<
                execution::is_execution_policy<ExecutionPolicy>()>>
  static void apply(ExecutionPolicy executionPolicy,
                    std::tuple<ArrayView...>&& arrayViews,
                    const Function& function) {
    apply(executionPolicy, std::move(arrayViews), detail::no_mask, function);
  }

  /// brief   Apply "For-Each" method
  ///
  /// details Apply ForEach with default execution policy and no mask.
  template <typename... ArrayView, typename Function>
  static void apply(std::tuple<ArrayView...>&& arrayViews,
                    const Function& function) {
    apply(execution::seq, std::move(arrayViews), function);
  }
};

/// brief   Construct ArrayForEach and call apply
///
/// details Construct an ArrayForEach<ItrDims...> using std::integer_sequence
///         <int, ItrDims...>. Remaining arguments are forwarded to apply
///         method.
template <int... ItrDims, typename... Args>
void arrayForEachDim(std::integer_sequence<int, ItrDims...>, Args&&... args) {
  ArrayForEach<ItrDims...>::apply(std::forward<Args>(args)...);
}

}  // namespace helpers
}  // namespace array
}  // namespace atlas
