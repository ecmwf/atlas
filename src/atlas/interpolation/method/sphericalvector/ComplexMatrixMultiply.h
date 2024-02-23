/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>
#include <iomanip>
#include <limits>
#include <sstream>
#include <tuple>
#include <type_traits>

#include "atlas/array/ArrayView.h"
#include "atlas/array/Range.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/interpolation/method/sphericalvector/Types.h"
#include "atlas/parallel/omp/omp.h"

// OpemMP support seems to be buggy for NVHPC
#ifdef __NVCOMPILER
#warning turning off OpenMP for nvhpc build
#undef atlas_omp_parallel_for
#define atlas_omp_parallel_for for
#endif

namespace atlas {
namespace interpolation {
namespace method {
namespace detail {

struct TwoVectorTag {};
struct ThreeVectorTag {};
constexpr auto twoVector = TwoVectorTag{};
constexpr auto threeVector = ThreeVectorTag{};

template <typename VectorTag>
constexpr bool isVectorTag() {
  return std::is_same_v<VectorTag, TwoVectorTag> ||
         std::is_same_v<VectorTag, ThreeVectorTag>;
}

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
  ComplexMatrixMultiply(const ComplexMatrix& complexWeights,
                        const RealMatrix& realWeights)
      : complexWeightsPtr_{&complexWeights}, realWeightsPtr_{&realWeights} {
    if constexpr (ATLAS_BUILD_TYPE_DEBUG) {
      ATLAS_ASSERT(complexWeightsPtr_->rows() == realWeightsPtr_->rows());
      ATLAS_ASSERT(complexWeightsPtr_->cols() == realWeightsPtr_->cols());
      ATLAS_ASSERT(complexWeightsPtr_->nonZeros() ==
                   realWeightsPtr_->nonZeros());
      for (auto rowIndex = Index{0}; rowIndex < complexWeightsPtr_->rows();
           ++rowIndex) {
        for (auto [complexRowIter, realRowIter] = rowIters(rowIndex);
             complexRowIter; ++complexRowIter, ++realRowIter) {
          ATLAS_ASSERT(realRowIter);
          ATLAS_ASSERT(realRowIter.row() == complexRowIter.row());
          ATLAS_ASSERT(realRowIter.col() == complexRowIter.col());
          // tinyNum ~= 2.3e-13 for double.
          constexpr auto tinyNum = 1024 * std::numeric_limits<Real>::epsilon();
          const auto complexMagnitude = std::abs(complexRowIter.value());
          const auto realValue = realRowIter.value();
          const auto error = std::abs(complexMagnitude - realValue);

          const auto printError = [&]() {
            auto strStream = std::ostringstream{};
            strStream << "Error complex components: ";
            strStream << std::setprecision(19) << error;
            return strStream.str();
          };
          ATLAS_ASSERT_MSG(error < tinyNum, printError());
        }
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
            typename = std::enable_if_t<isVectorTag<VectorType>()>>
  void apply(const array::ArrayView<const Value, Rank>& sourceView,
             array::ArrayView<Value, Rank>& targetView, VectorType) const {
    if constexpr (std::is_same_v<VectorType, TwoVectorTag>) {
      applyTwoVector(sourceView, targetView);
    }

    if constexpr (std::is_same_v<VectorType, ThreeVectorTag>) {
      applyThreeVector(sourceView, targetView);
    }

    return;
  }

 private:
  /// @brief Apply 2-vector MatMul.
  template <typename Value, int Rank>
  void applyTwoVector(const array::ArrayView<const Value, Rank>& sourceView,
                      array::ArrayView<Value, Rank>& targetView) const {
    // We could probably optimise contiguous arrays using
    // reinterpret_cast<std::complex<double>*>(view.data()). This is fine
    // according to the C++ standard!
    atlas_omp_parallel_for(auto rowIndex = Index{0};
                           rowIndex < complexWeightsPtr_->rows(); ++rowIndex) {
      auto targetSlice = sliceColumn(targetView, rowIndex);
      if constexpr (InitialiseTarget) {
        targetSlice.assign(0.);
      }

      for (auto complexRowIter = complexWeightsPtr_->rowIter(rowIndex);
           complexRowIter; ++complexRowIter) {
        const auto colIndex = complexRowIter.col();
        const auto complexWeight = complexRowIter.value();
        const auto sourceSlice = sliceColumn(sourceView, colIndex);

        array::helpers::arrayForEachDim(
            slicedColumnDims<Rank>{}, std::tie(sourceSlice, targetSlice),
            [&](auto&& sourceElem, auto&& targetElem) {
              const auto targetVector =
                  complexWeight * Complex(sourceElem(0), sourceElem(1));
              targetElem(0) += targetVector.real();
              targetElem(1) += targetVector.imag();
            });
      }
    }
  }

  /// @brief Apply 3-vector MatMul.
  template <typename Value, int Rank>
  void applyThreeVector(const array::ArrayView<const Value, Rank>& sourceView,
                        array::ArrayView<Value, Rank>& targetView) const {
    atlas_omp_parallel_for(auto rowIndex = Index{0};
                           rowIndex < complexWeightsPtr_->rows(); ++rowIndex) {
      auto targetSlice = sliceColumn(targetView, rowIndex);
      if constexpr (InitialiseTarget) {
        targetSlice.assign(0.);
      }

      for (auto [complexRowIter, realRowIter] = rowIters(rowIndex);
           complexRowIter; ++complexRowIter, ++realRowIter) {
        const auto colIndex = complexRowIter.col();
        const auto complexWeight = complexRowIter.value();
        const auto realWeight = realRowIter.value();
        const auto sourceSlice = sliceColumn(sourceView, colIndex);

        array::helpers::arrayForEachDim(
            slicedColumnDims<Rank>{}, std::tie(sourceSlice, targetSlice),
            [&](auto&& sourceElem, auto&& targetElem) {
              const auto targetVector =
                  complexWeight * Complex(sourceElem(0), sourceElem(1));
              targetElem(0) += targetVector.real();
              targetElem(1) += targetVector.imag();
              targetElem(2) += realWeight * sourceElem(2);
            });
      }
    }
  }

  /// @brief Return a pair of complex and real row iterators
  std::pair<ComplexMatrix::RowIter, RealMatrix::RowIter> rowIters(
      Index rowIndex) const {
    return std::make_pair(complexWeightsPtr_->rowIter(rowIndex),
                          realWeightsPtr_->rowIter(rowIndex));
  }

  /// @brief Makes the slice arrayView.slice(index, Range::all()...).
  template <typename View, typename Index>
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

  const ComplexMatrix* complexWeightsPtr_{};
  const RealMatrix* realWeightsPtr_{};
};

}  // namespace detail
}  // namespace method
}  // namespace interpolation
}  // namespace atlas
