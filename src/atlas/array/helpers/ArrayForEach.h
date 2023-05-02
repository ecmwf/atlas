/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>
#include <tuple>
#include <type_traits>

#include "atlas/array/ArrayView.h"
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
std::string policy_name() {
  return "unsupported";
};
template <>
std::string policy_name<sequenced_policy>() {
  return "sequenced_policy";
};
template <>
std::string policy_name<unsequenced_policy>() {
  return "unsequenced_policy";
};
template <>
std::string policy_name<parallel_unsequenced_policy>() {
  return "parallel_unsequenced_policy";
};

template <typename execution_policy>
std::string policy_name(execution_policy) {
  return policy_name<execution_policy>();
}

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

template <typename... SlicerArgs, template <typename, int> typename View,
          typename Value, int Rank, typename... ArrayViews>
auto makeSlices(const std::tuple<SlicerArgs...>& slicerArgs,
                View<Value, Rank>& arrayView, ArrayViews&... arrayViews) {

  // "Lambdafy" slicer apply method to work with std::apply.
  const auto slicer = [&arrayView](const auto&... args) {
    return std::make_tuple(arrayView.slice(args...));
  };

  // Fill out the remaining slicerArgs with Range::all().
  constexpr auto Dim = sizeof...(SlicerArgs);
  const auto paddedArgs = std::tuple_cat(slicerArgs, argPadding<Rank - Dim>());

  if constexpr(sizeof...(ArrayViews) > 0) {

    // Recurse until all views are sliced.
    return std::tuple_cat(std::apply(slicer, paddedArgs),
                          makeSlices(slicerArgs, arrayViews...));
  }
  else {
    return std::apply(slicer, paddedArgs);
  }
}

template <typename ExecutionPolicy, int Dim, int... ItrDims>
struct ArrayForEachImpl;

template <typename ExecutionPolicy, int Dim, int ItrDim, int... ItrDims>
struct ArrayForEachImpl<ExecutionPolicy, Dim, ItrDim, ItrDims...> {
  template <typename... ArrayViews, typename Mask, typename Function,
            typename... SlicerArgs, typename... MaskArgs>
  static void apply(std::tuple<ArrayViews...>& arrayViews,
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
              arrayViews, mask, function,
              tuplePushBack(slicerArgs, idx),
              tuplePushBack(maskArgs, idx));
        }
        else {
          // Retain current execution policy.
          ArrayForEachImpl<ExecutionPolicy, Dim + 1, ItrDims...>::apply(
              arrayViews, mask, function,
              tuplePushBack(slicerArgs, idx),
              tuplePushBack(maskArgs, idx));
        }
      });
    }
    // Add a RangeAll to arguments.
    else {
      ArrayForEachImpl<ExecutionPolicy, Dim + 1, ItrDim, ItrDims...>::apply(
          arrayViews, mask, function,
          tuplePushBack(slicerArgs, Range::all()),
          maskArgs);
    }
  }
};

template <typename ExecutionPolicy, int Dim>
struct ArrayForEachImpl<ExecutionPolicy, Dim> {
  template <typename... ArrayViews, typename Mask, typename Function,
            typename... SlicerArgs, typename... MaskArgs>
  static void apply(std::tuple<ArrayViews...>& arrayViews,
                    const Mask& mask,
                    const Function& function,
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
  template <typename... ArrayViews, typename Mask, typename Function>
  static void apply(const eckit::Parametrisation& conf,
                    const std::tuple<ArrayViews...>& arrayViews,
                    const Mask& mask, const Function& function) {

    auto execute = [&](auto execution_policy) {
        // Make a copy of views to simplify constness and forwarding.
        auto arrayViewsCopy = arrayViews;

        detail::ArrayForEachImpl<decltype(execution_policy), 0, ItrDims...>::apply(
            arrayViewsCopy, mask, function, std::make_tuple(), std::make_tuple());
    };

    using namespace execution;
    std::string execution_policy;
    if( conf.get("execution_policy",execution_policy) ) {
        if (execution_policy == policy_name(par_unseq)) {
            execute(par_unseq);
        }
        else if (execution_policy == policy_name(par)) {
            execute(par);
        }
        else if (execution_policy == policy_name(unseq)) {
            execute(unseq);
        }
        else if (execution_policy == policy_name(seq)) {
            execute(seq);
        }
    }
    else {
      execute(par_unseq);
    }
  }

  template <typename ExecutionPolicy, typename... ArrayViews, typename Mask,
            typename Function, std::enable_if_t<!std::is_base_of_v<eckit::Parametrisation,ExecutionPolicy>,int> =0>
  static void apply(ExecutionPolicy, const std::tuple<ArrayViews...>& arrayViews, const Mask& mask, const Function& function) {
      apply(option::execution_policy<ExecutionPolicy>(), arrayViews, mask, function);
  }

  template <typename... ArrayViews, typename Mask, typename Function>
  static void apply(const std::tuple<ArrayViews...>& arrayViews, const Mask& mask, const Function& function) {
      apply(util::NoConfig(), arrayViews, mask, function);
  }

  /// brief   Apply "For-Each" method.
  ///
  /// details Apply "For-Each" without a mask.
  template <typename... ArrayViews, typename Function>
  static void apply(const eckit::Parametrisation& conf, const std::tuple<ArrayViews...>& arrayViews, const Function& function) {
    constexpr auto no_mask = [](auto args...) { return 0; };
    apply(conf, arrayViews, no_mask, function);
  }

  template <typename ExecutionPolicy, typename... ArrayViews,
            typename Function, std::enable_if_t<!std::is_base_of_v<eckit::Parametrisation,ExecutionPolicy>,int> =0>
  static void apply(ExecutionPolicy, const std::tuple<ArrayViews...>& arrayViews, const Function& function) {
    apply(option::execution_policy<ExecutionPolicy>(), arrayViews, function);
  }

  template <typename... ArrayViews, typename Function>
  static void apply(const std::tuple<ArrayViews...>& arrayViews, const Function& function) {
    apply(util::NoConfig{}, arrayViews, function);
  }

};

}  // namespace helpers
}  // namespace array
}  // namespace atlas
