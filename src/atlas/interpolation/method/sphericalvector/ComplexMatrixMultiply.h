#pragma once

#include <array>
#include <memory>
#include <tuple>
#include <type_traits>

#include "atlas/array/ArrayView.h"
#include "atlas/array/Range.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/interpolation/method/sphericalvector/SphericalVector.h"
#include "atlas/parallel/omp/omp.h"

namespace atlas {
namespace interpolation {
namespace method {
namespace detail {

using Complex = SphericalVector::Complex;
using ComplexMatPtr = SphericalVector::ComplexMatPtr;
using RealMatPtr = SphericalVector::RealMatPtr;

struct TwoVectorTag {};
struct ThreeVectorTag {};
constexpr TwoVectorTag twoVector{};
constexpr ThreeVectorTag threeVector{};

template <typename VectorTag>
using IsVectorTag = std::enable_if_t<std::is_same_v<VectorTag, TwoVectorTag> ||
                                     std::is_same_v<VectorTag, ThreeVectorTag>>;

template <bool InitialiseTarget>
class ComplexMatrixMultiply {
 public:
  ComplexMatrixMultiply() = default;
  ComplexMatrixMultiply(const ComplexMatPtr& complexWeights,
                        const RealMatPtr& realWeights);

  template <typename Value, int Rank, typename VectorType,
            typename = IsVectorTag<VectorType>>
  void apply(const array::ArrayView<const Value, Rank>& sourceView,
             array::ArrayView<Value, Rank>& targetView, VectorType) const;

 private:
  template <typename View, typename Index>
  static auto sliceColumn(View& arrayView, Index index);

  template <int Rank>
  using slicedColumnDims = std::make_integer_sequence<int, Rank - 2>;

  const ComplexMatPtr complexWeights_{};
  const RealMatPtr realWeights_{};
};

template <bool InitialiseTarget>
inline ComplexMatrixMultiply<InitialiseTarget>::ComplexMatrixMultiply(
    const ComplexMatPtr& complexWeights, const RealMatPtr& realWeights)
    : complexWeights_{complexWeights}, realWeights_{realWeights} {}

template <bool InitialiseTarget>
template <typename Value, int Rank, typename VectorType, typename>
void ComplexMatrixMultiply<InitialiseTarget>::apply(
    const array::ArrayView<const Value, Rank>& sourceView,
    array::ArrayView<Value, Rank>& targetView, VectorType) const {

  const auto* outerIndices = complexWeights_->outerIndexPtr();
  const auto* innerIndices = complexWeights_->innerIndexPtr();
  const auto* complexWeightValues = complexWeights_->valuePtr();
  const auto* realWeightValues = realWeights_->valuePtr();
  const auto nRows = complexWeights_->outerSize();

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

template <bool InitialiseTarget>
template <typename View, typename Index>
auto ComplexMatrixMultiply<InitialiseTarget>::sliceColumn(View& arrayView,
                                                          Index index) {

  constexpr auto Rank = std::decay_t<View>::rank();
  using RangeAll = decltype(array::Range::all());

  const auto slicerArgs =
      std::tuple_cat(std::make_tuple(index), std::array<RangeAll, Rank - 1>{});
  const auto slicer = [&](auto&&... args) { return arrayView.slice(args...); };

  return std::apply(slicer, slicerArgs);
}

}  // namespace detail
}  // namespace method
}  // namespace interpolation
}  // namespace atlas
