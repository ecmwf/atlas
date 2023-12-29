#pragma once

#include <array>
#include <tuple>

#include "atlas/array/ArrayView.h"
#include "atlas/array/Range.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/interpolation/method/sphericalvector/SphericalVector.h"

namespace atlas {
namespace interpolation {
namespace method {
namespace detail {

using Complex = SphericalVector::Complex;
using ComplexMatrix = SphericalVector::ComplexMatrix;
using RealMatrix = SphericalVector::RealMatrix;

class ComplexMatrixMultiply {
 public:
  ComplexMatrixMultiply(const ComplexMatrix& complexWeights,
                        const RealMatrix& realWeights);

  template <bool InitialiseTarget = true, typename Value, int Rank>
  void ApplyTwoVector(const array::ArrayView<const Value, Rank>& sourceView,
                      array::ArrayView<Value, Rank>& targetView) const;

 private:
  template <int Rank>
  using columnDims = std::make_integer_sequence<int, Rank - 2>;

  template <typename View, typename Index>
  static auto sliceColumn(View& arrayView, Index index);

  const ComplexMatrix& complexWeights_;
  const RealMatrix& realWeights_;
};

inline ComplexMatrixMultiply::ComplexMatrixMultiply(
    const ComplexMatrix& complexWeights, const RealMatrix& realWeights)
    : complexWeights_{complexWeights}, realWeights_{realWeights} {}

template <bool InitialiseTarget, typename Value, int Rank>
void ComplexMatrixMultiply::ApplyTwoVector(
    const array::ArrayView<const Value, Rank>& sourceView,
    array::ArrayView<Value, Rank>& targetView) const {
  using Index = ComplexMatrix::Index;
  using InnerIterator = ComplexMatrix::InnerIterator;

  for (auto k = Index{}; k < complexWeights_.outerSize(); ++k) {
    auto it = InnerIterator(complexWeights_, k);
    const auto rowIndex = it.row();
    auto targetSlice = sliceColumn(targetView, rowIndex);
    if constexpr (InitialiseTarget) {
      targetSlice.assign(Value{});
    }

    do {
      const auto colIndex = it.col();
      const auto& weight = it.value();
      const auto sourceSlice = sliceColumn(sourceView, colIndex);

      array::helpers::arrayForEachDim(
          columnDims<Rank>{}, std::tie(sourceSlice, targetSlice),
          [&](auto&& sourceElem, auto&& targetElem) {
            const auto targetVector =
                weight * Complex(sourceElem(0), sourceElem(1));
            targetElem(0) += targetVector.real();
            targetElem(1) += targetVector.imag();
          });

    } while (++it);
  }
}

template <typename View, typename Index>
auto ComplexMatrixMultiply::sliceColumn(View& arrayView, Index index) {
  constexpr auto Rank = std::decay_t<View>::rank();
  using RangeAll = decltype(array::Range::all());
  const auto slicerArgs =
      std::tuple_cat(std::make_tuple(index), std::array<RangeAll, Rank - 1>{});
  return std::apply([&](auto&&... args) { return arrayView.slice(args...); },
                    slicerArgs);
}

}  // namespace detail
}  // namespace method
}  // namespace interpolation
}  // namespace atlas
