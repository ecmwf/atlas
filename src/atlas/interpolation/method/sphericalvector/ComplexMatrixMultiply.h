/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>
#include <memory>
#include <tuple>
#include <type_traits>

#include "atlas/array/ArrayView.h"
#include "atlas/array/Range.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/interpolation/method/sphericalvector/Types.h"
#include "atlas/parallel/omp/omp.h"

namespace atlas {
namespace interpolation {
namespace method {
namespace detail {

using ComplexMatPtr = std::shared_ptr<ComplexMatrix>;
using RealMatPtr = std::shared_ptr<RealMatrix>;

struct TwoVectorTag {};
struct ThreeVectorTag {};
constexpr TwoVectorTag twoVector{};
constexpr ThreeVectorTag threeVector{};

template <typename VectorTag>
using IsVectorTag = std::enable_if_t<std::is_same_v<VectorTag, TwoVectorTag> ||
                                     std::is_same_v<VectorTag, ThreeVectorTag>>;

/// @brief   Helper class to perform complex matrix multiplications
///
/// @details Performs matrix multiplication between fields of 2-vectors or
///          3-vectors. Fields must have Rank >= 2. Here, the assumption is
///          that Dim = 0 is the horizontal dimension, and Dim = (Rank - 1) is
///          the vector element dimension.
template <bool InitialiseTarget>
class ComplexMatrixMultiply {
 public:
  ComplexMatrixMultiply() = default;

  /// @brief   Construct object from sparse matrices.
  ///
  /// @details complexWeights is a SparseMatrix of weights. realWeights is a
  ///          SparseMatrix containing the magnitudes of the elements of
  ///          complexWeights.
  ComplexMatrixMultiply(const ComplexMatPtr& complexWeights,
                        const RealMatPtr& realWeights)
      : complexWeights_{complexWeights}, realWeights_{realWeights} {

    if constexpr (ATLAS_BUILD_TYPE_DEBUG) {
      ATLAS_ASSERT(complexWeights->rows() == realWeights->rows());
      ATLAS_ASSERT(complexWeights->cols() == realWeights->cols());
      ATLAS_ASSERT(complexWeights->nonZeros() == realWeights->nonZeros());

      for (auto i = Size{0}; i < complexWeights->rows() + 1; ++i) {
        ATLAS_ASSERT(complexWeights_->outer()[i] == realWeights_->outer()[i]);
      }

      for (auto i = Size{0}; i < complexWeights->nonZeros(); ++i) {
        ATLAS_ASSERT(complexWeights_->inner()[i] == realWeights_->inner()[i]);
      }
    }
  }

  /// @brief   Apply complex matrix vector multiplication.
  ///
  /// @details Multiply weights by the elements in sourceView to give
  ///          elements in targetView. If VectorType == TwoVectorTag,
  ///          complexWeights are applied to the horizontal elements of
  ///          sourceView. If VectorType == ThreeVectorTag, then realWeights
  ///          are additionally applied to the vertical elements of sourceView.
  template <typename Value, int Rank, typename VectorType,
            typename = IsVectorTag<VectorType>>
  void apply(const array::ArrayView<const Value, Rank>& sourceView,
             array::ArrayView<Value, Rank>& targetView, VectorType) const {

    const auto* outerIndices = complexWeights_->outer();
    const auto* innerIndices = complexWeights_->inner();
    const auto* complexWeightValues = complexWeights_->data();
    const auto* realWeightValues = realWeights_->data();
    const auto nRows = complexWeights_->rows();

    using Index = std::decay_t<decltype(*innerIndices)>;

    atlas_omp_parallel_for(auto rowIndex = Index{0}; rowIndex < nRows;
                           ++rowIndex) {

      auto targetSlice = sliceColumn(targetView, rowIndex);
      if constexpr (InitialiseTarget) { targetSlice.assign(0.); }
      for (auto dataIndex = outerIndices[rowIndex];
           dataIndex < outerIndices[rowIndex + 1]; ++dataIndex) {

        const auto colIndex = innerIndices[dataIndex];
        const auto sourceSlice = sliceColumn(sourceView, colIndex);

        array::helpers::arrayForEachDim(
            slicedColumnDims<Rank>{}, std::tie(sourceSlice, targetSlice),
            [&](auto&& sourceElem, auto&& targetElem) {
              const auto targetVector = complexWeightValues[dataIndex] *
                                        Complex(sourceElem(0), sourceElem(1));
              targetElem(0) += targetVector.real();
              targetElem(1) += targetVector.imag();

              if constexpr (std::is_same_v<VectorType, ThreeVectorTag>) {
                  targetElem(2) += realWeightValues[dataIndex] * sourceElem(2);
              }
            });
      }
    }
  }

 private:
  template <typename View, typename Index>

  /// @brief Makes the slice arrayView.slice(index, Range::all()...).
  static auto sliceColumn(View& arrayView, Index index) {
    constexpr auto Rank = std::decay_t<View>::rank();
    using RangeAll = decltype(array::Range::all());

    const auto slicerArgs = std::tuple_cat(std::make_tuple(index),
                                           std::array<RangeAll, Rank - 1>{});
    const auto slicer = [&](auto&&... args) {
      return arrayView.slice(args...);
    };

    return std::apply(slicer, slicerArgs);
  }

  /// @brief Defines a sequence of iteration Dims for a sliced column.
  template <int Rank>
  using slicedColumnDims = std::make_integer_sequence<int, Rank - 2>;

  ComplexMatPtr complexWeights_{};
  RealMatPtr realWeights_{};
};

}  // namespace detail
}  // namespace method
}  // namespace interpolation
}  // namespace atlas
