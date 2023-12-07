/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/library/defines.h"

#include <cmath>
#include <tuple>
#include <utility>

#include "atlas/array/ArrayView.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/array/Range.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/interpolation/Cache.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/interpolation/method/sphericalvector/SphericalVector.h"
#include "atlas/linalg/sparse.h"
#include "atlas/option/Options.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Constants.h"
#include "atlas/util/UnitSphere.h"

#include "eckit/linalg/Triplet.h"

namespace atlas {
namespace interpolation {
namespace method {

using Complex = SphericalVector::Complex;

#if ATLAS_HAVE_EIGEN
template <typename Value>
using SparseMatrix = SphericalVector::SparseMatrix<Value>;
using RealMatrixMap = Eigen::Map<const SparseMatrix<double>>;
using ComplexTriplets = std::vector<Eigen::Triplet<Complex>>;
using RealTriplets = std::vector<Eigen::Triplet<double>>;
#endif

using EckitMatrix = eckit::linalg::SparseMatrix;

namespace {

MethodBuilder<SphericalVector> __builder("spherical-vector");

#if ATLAS_HAVE_EIGEN
RealMatrixMap makeMatrixMap(const EckitMatrix& baseMatrix) {
  return RealMatrixMap(baseMatrix.rows(), baseMatrix.cols(),
                       baseMatrix.nonZeros(), baseMatrix.outer(),
                       baseMatrix.inner(), baseMatrix.data());
}

template <typename Matrix, typename Index>
auto getInnerIt(const Matrix& matrix, Index k) {
  return typename Matrix::InnerIterator(matrix, k);
}

template <typename Functor, typename Matrix>
void sparseMatrixForEach(const Functor& functor, const Matrix& matrix) {

  using Index = decltype (matrix.outerSize());
  atlas_omp_parallel_for (auto k = Index{}; k < matrix.outerSize(); ++k) {
    for (auto it = getInnerIt(matrix, k); it; ++it) {
      functor(it.row(), it.col(), it.value());
    }
  }
}

template <typename Functor, typename Matrix1, typename Matrix2>
void sparseMatrixForEach(const Functor& functor, const Matrix1& matrix1,
                         const Matrix2& matrix2) {
  using Index = decltype (matrix1.outerSize());
  atlas_omp_parallel_for (auto k = Index{}; k < matrix1.outerSize(); ++k) {
    for (auto [it1, it2] =
             std::make_pair(getInnerIt(matrix1, k), getInnerIt(matrix2, k));
         it1; ++it1, ++it2) {
      functor(it1.row(), it1.col(), it1.value(), it2.value());
    }
  }
}

template <typename SourceView, typename TargetView, typename Functor,
          typename... Matrices>
void matrixMultiply(const SourceView& sourceView, TargetView& targetView,
                    const Functor& multiplyFunctor,
                    const Matrices&... matrices) {

  const auto multiplyColumn = [&](auto i, auto j, const auto&... weights) {
    constexpr auto Rank = std::decay_t<SourceView>::rank();
    if constexpr (Rank == 2) {
        const auto sourceSlice = sourceView.slice(j, array::Range::all());
        auto targetSlice = targetView.slice(i, array::Range::all());
        multiplyFunctor(sourceSlice, targetSlice, weights...);
    } else if constexpr(Rank == 3) {
        const auto sourceSlice =
            sourceView.slice(j, array::Range::all(), array::Range::all());
        auto targetSlice =
            targetView.slice(i, array::Range::all(), array::Range::all());
        array::helpers::ArrayForEach<0>::apply(
            std::tie(sourceSlice, targetSlice),
            [&](auto&& sourceVars, auto&& targetVars) {
              multiplyFunctor(sourceVars, targetVars, weights...);
            });
    } else {
      ATLAS_NOTIMPLEMENTED;
    }
  };

  sparseMatrixForEach(multiplyColumn, matrices...);
}
#endif

}  // namespace

void SphericalVector::do_setup(const Grid& source, const Grid& target,
                               const Cache&) {
  ATLAS_NOTIMPLEMENTED;
}

void SphericalVector::do_setup(const FunctionSpace& source,
                               const FunctionSpace& target) {
  ATLAS_TRACE("interpolation::method::SphericalVector::do_setup");
  source_ = source;
  target_ = target;

  if (target_.size() == 0) {
    return;
  }

#if ATLAS_HAVE_EIGEN

  setMatrix(Interpolation(interpolationScheme_, source_, target_));

  // Get matrix data.
  const auto nRows = matrix().rows();
  const auto nCols = matrix().cols();
  const auto nNonZeros = matrix().nonZeros();
  const auto baseWeights = makeMatrixMap(matrix());

  // Note: need to store copy of weights as Eigen3 sorts compressed rows by j
  // whereas eckit does not.
  complexWeights_ = std::make_shared<ComplexMatrix>(nRows, nCols);
  realWeights_ = std::make_shared<RealMatrix>(nRows, nCols);
  auto complexTriplets = ComplexTriplets(nNonZeros);
  auto realTriplets = RealTriplets(nNonZeros);

  const auto sourceLonLats = array::make_view<double, 2>(source_.lonlat());
  const auto targetLonLats = array::make_view<double, 2>(target_.lonlat());

  const auto setComplexWeights = [&](auto i, auto j, const auto& weight) {
    const auto sourceLonLat =
        PointLonLat(sourceLonLats(j, 0), sourceLonLats(j, 1));
    const auto targetLonLat =
        PointLonLat(targetLonLats(i, 0), targetLonLats(i, 1));

    const auto alpha = util::greatCircleCourse(sourceLonLat, targetLonLat);

    const auto deltaAlpha =
        (alpha.first - alpha.second) * util::Constants::degreesToRadians();

    const auto idx = std::distance(baseWeights.valuePtr(), &weight);

    complexTriplets[idx] = {int(i), int(j), std::polar(weight, deltaAlpha)};
    realTriplets[idx] = {int(i), int(j), weight};
  };

  sparseMatrixForEach(setComplexWeights, baseWeights);
  complexWeights_->setFromTriplets(complexTriplets.begin(),
                                   complexTriplets.end());
  realWeights_->setFromTriplets(realTriplets.begin(), realTriplets.end());

  ATLAS_ASSERT(complexWeights_->nonZeros() == matrix().nonZeros());

#else
  ATLAS_THROW_EXCEPTION("atlas has been compiled without Eigen");
#endif
}

void SphericalVector::print(std::ostream&) const { ATLAS_NOTIMPLEMENTED; }

void SphericalVector::do_execute(const FieldSet& sourceFieldSet,
                                 FieldSet& targetFieldSet,
                                 Metadata& metadata) const {
  ATLAS_TRACE("atlas::interpolation::method::SphericalVector::do_execute()");

  const auto nFields = sourceFieldSet.size();
  ATLAS_ASSERT(nFields == targetFieldSet.size());

  for (auto i = 0; i < sourceFieldSet.size(); ++i) {
    do_execute(sourceFieldSet[i], targetFieldSet[i], metadata);
  }
}

void SphericalVector::do_execute(const Field& sourceField, Field& targetField,
                                 Metadata&) const {
  ATLAS_TRACE("atlas::interpolation::method::SphericalVector::do_execute()");

  const auto fieldType = sourceField.metadata().getString("type", "");
  if (fieldType != "vector") {

    auto metadata = Metadata();
    Method::do_execute(sourceField, targetField, metadata);

    return;
  }

  if (target_.size() == 0) {
    return;
  }

  ATLAS_ASSERT_MSG(sourceField.variables() == 2 || sourceField.variables() == 3,
                   "Vector field can only have 2 or 3 components.");

  Method::check_compatibility(sourceField, targetField, matrix());

  haloExchange(sourceField);

  if (sourceField.datatype().kind() == array::DataType::KIND_REAL64) {
    interpolate_vector_field<double>(sourceField, targetField);
  } else if (sourceField.datatype().kind() == array::DataType::KIND_REAL32) {
    interpolate_vector_field<float>(sourceField, targetField);
  } else {
    ATLAS_NOTIMPLEMENTED;
  }

  targetField.set_dirty();
}

template <typename Value>
void SphericalVector::interpolate_vector_field(const Field& sourceField,
                                               Field& targetField) const {
  if (sourceField.rank() == 2) {
    interpolate_vector_field<Value, 2>(sourceField, targetField);
  } else if (sourceField.rank() == 3) {
    interpolate_vector_field<Value, 3>(sourceField, targetField);
  } else {
    ATLAS_NOTIMPLEMENTED;
  }
}

template <typename Value, int Rank>
void SphericalVector::interpolate_vector_field(const Field& sourceField,
                                               Field& targetField) const {

  const auto sourceView = array::make_view<Value, Rank>(sourceField);
  auto targetView = array::make_view<Value, Rank>(targetField);
  targetView.assign(0.);

#if ATLAS_HAVE_EIGEN
  const auto horizontalComponent = [](const auto& sourceVars, auto& targetVars,
                                      const auto& complexWeight) {
    const auto sourceVector = Complex(sourceVars(0), sourceVars(1));
    const auto targetVector = complexWeight * sourceVector;
    targetVars(0) += targetVector.real();
    targetVars(1) += targetVector.imag();
  };

  if (sourceField.variables() == 2) {
    matrixMultiply(sourceView, targetView, horizontalComponent,
                   *complexWeights_);
    return;
  } else if (sourceField.variables() == 3) {

    const auto horizontalAndVerticalComponent = [&](
        const auto& sourceVars, auto& targetVars, const auto& complexWeight,
        const auto& realWeight) {
      horizontalComponent(sourceVars, targetVars, complexWeight);
      targetVars(2) += realWeight * sourceVars(2);
    };

    matrixMultiply(sourceView, targetView, horizontalAndVerticalComponent,
                   *complexWeights_, *realWeights_);

    return;
  }
#else
  ATLAS_THROW_EXCEPTION("atlas has been compiled without Eigen");
#endif

  ATLAS_NOTIMPLEMENTED;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
