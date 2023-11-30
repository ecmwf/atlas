/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/library/defines.h"
#if ATLAS_HAVE_EIGEN

#include <cmath>
#include <tuple>

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

template <typename Value>
using SparseMatrix = SphericalVector::SparseMatrix<Value>;
using RealMatrixMap = Eigen::Map<const SparseMatrix<double>>;
using ComplexTriplets = std::vector<Eigen::Triplet<Complex>>;
using EckitMatrix = eckit::linalg::SparseMatrix;

namespace {

MethodBuilder<SphericalVector> __builder("spherical-vector");

RealMatrixMap makeMatrixMap(const EckitMatrix& baseMatrix) {
  return RealMatrixMap(baseMatrix.rows(), baseMatrix.cols(),
                       baseMatrix.nonZeros(), baseMatrix.outer(),
                       baseMatrix.inner(), baseMatrix.data());
}

template <typename MatrixT, typename Functor>
void sparseMatrixForEach(const MatrixT& matrix, const Functor& functor) {

  atlas_omp_parallel_for (auto k = 0; k < matrix.outerSize(); ++k) {
    for (auto it = typename MatrixT::InnerIterator(matrix, k); it; ++it) {
          functor(it.value(), it.row(), it.col());
    }
  }
}

template <typename MatrixT, typename SourceView, typename TargetView,
          typename Functor>
void matrixMultiply(const MatrixT& matrix, SourceView&& sourceView,
                    TargetView&& targetView, const Functor& mappingFunctor) {

  sparseMatrixForEach(matrix, [&](const auto& weight, auto i, auto j) {

    constexpr auto rank = std::decay_t<SourceView>::rank();
    if constexpr(rank == 2) {
        const auto sourceSlice = sourceView.slice(j, array::Range::all());
        auto targetSlice = targetView.slice(i, array::Range::all());
        mappingFunctor(weight, sourceSlice, targetSlice);
    }
    else if constexpr(rank == 3) {
        const auto iterationFuctor = [&](auto&& sourceVars, auto&& targetVars) {
          mappingFunctor(weight, sourceVars, targetVars);
        };
        const auto sourceSlice =
            sourceView.slice(j, array::Range::all(), array::Range::all());
        auto targetSlice =
            targetView.slice(i, array::Range::all(), array::Range::all());
        array::helpers::ArrayForEach<0>::apply(
            std::tie(sourceSlice, targetSlice), iterationFuctor);
    }
    else {
      ATLAS_NOTIMPLEMENTED;
    }
  });
}

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

  const auto baseInterpolator =
      Interpolation(interpolationScheme_, source_, target_);
  setMatrix(MatrixCache(baseInterpolator));

  // Get matrix data.
  const auto nRows = matrix().rows();
  const auto nCols = matrix().cols();
  const auto nNonZeros = matrix().nonZeros();
  const auto realWeights = makeMatrixMap(matrix());

  complexWeights_ = std::make_shared<ComplexMatrix>(nRows, nCols);
  auto complexTriplets = ComplexTriplets(nNonZeros);

  const auto sourceLonLats = array::make_view<double, 2>(source_.lonlat());
  const auto targetLonLats = array::make_view<double, 2>(target_.lonlat());

  sparseMatrixForEach(realWeights, [&](const auto& weight, auto i, auto j) {

    const auto sourceLonLat =
        PointLonLat(sourceLonLats(j, 0), sourceLonLats(j, 1));
    const auto targetLonLat =
        PointLonLat(targetLonLats(i, 0), targetLonLats(i, 1));

    const auto alpha = util::greatCircleCourse(sourceLonLat, targetLonLat);

    const auto deltaAlpha =
        (alpha.first - alpha.second) * util::Constants::degreesToRadians();

    const auto idx = std::distance(realWeights.valuePtr(), &weight);

    complexTriplets[idx] = {int(i), int(j), std::polar(weight, deltaAlpha)};
  });
  complexWeights_->setFromTriplets(complexTriplets.begin(),
                                   complexTriplets.end());

  ATLAS_ASSERT(complexWeights_->nonZeros() == matrix().nonZeros());
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

  const auto horizontalComponent = [](const auto& weight, auto&& sourceVars,
                                      auto&& targetVars) {
    const auto sourceVector = Complex(sourceVars(0), sourceVars(1));
    const auto targetVector = weight * sourceVector;
    targetVars(0) += targetVector.real();
    targetVars(1) += targetVector.imag();
  };

  matrixMultiply(*complexWeights_, sourceView, targetView, horizontalComponent);

  if (sourceField.variables() == 2) return;

  const auto verticalComponent = [](
      const auto& weight, auto&& sourceVars,
      auto&& targetVars) { targetVars(2) += weight * sourceVars(2); };

  const auto realWeights = makeMatrixMap(matrix());
  matrixMultiply(realWeights, sourceView, targetView, verticalComponent);
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas

#endif
