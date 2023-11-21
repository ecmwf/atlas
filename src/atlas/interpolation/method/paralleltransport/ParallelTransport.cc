/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

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
#include "atlas/interpolation/method/paralleltransport/ParallelTransport.h"
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

namespace {
MethodBuilder<ParallelTransport> __builder("parallel-transport");

template <typename MatrixT, typename Functor>
void spaceMatrixForEach(MatrixT&& matrix, const Functor& functor) {

  const auto nRows = matrix.rows();
  const auto nCols = matrix.cols();
  const auto rowIndices = matrix.outer();
  const auto colIndices = matrix.inner();
  auto valData = matrix.data();

  atlas_omp_parallel_for(auto i = size_t{}; i < nRows; ++i) {
    for (auto dataIdx = rowIndices[i]; dataIdx < rowIndices[i + 1]; ++dataIdx) {
      const auto j = size_t(colIndices[dataIdx]);
      auto&& value = valData[dataIdx];

      if
        constexpr(
            std::is_invocable_v<Functor, decltype(value), size_t, size_t>) {
          functor(value, i, j);
        }
      else if
        constexpr(std::is_invocable_v<Functor, decltype(value), size_t, size_t,
                                      size_t>) {
          functor(value, i, j, dataIdx);
        }
      else {
        ATLAS_NOTIMPLEMENTED;
      }
    }
  }
}

template <typename MatrixT, typename SourceView, typename TargetView,
          typename Functor>
void matrixMultiply(const MatrixT& matrix, SourceView&& sourceView,
                    TargetView&& targetView, const Functor& mappingFunctor) {

  spaceMatrixForEach(matrix, [&](const auto& weight, auto i, auto j) {

    constexpr auto rank = std::decay_t<decltype(sourceView)>::rank();
    if
      constexpr(rank == 2) {
        const auto sourceSlice = sourceView.slice(j, array::Range::all());
        auto targetSlice = targetView.slice(i, array::Range::all());
        mappingFunctor(weight, sourceSlice, targetSlice);
      }
    else if
      constexpr(rank == 3) {
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

void ParallelTransport::do_setup(const Grid& source, const Grid& target,
                                 const Cache&) {
  ATLAS_NOTIMPLEMENTED;
}

void ParallelTransport::do_setup(const FunctionSpace& source,
                                 const FunctionSpace& target) {
  ATLAS_TRACE("interpolation::method::ParallelTransport::do_setup");
  source_ = source;
  target_ = target;

  if (target_.size() == 0) {
    return;
  }

  const auto baseInterpolator =
      Interpolation(interpolationScheme_, source_, target_);
  setMatrix(MatrixCache(baseInterpolator));

  // Get matrix dimensions.
  const auto nRows = matrix().rows();
  const auto nCols = matrix().cols();
  const auto nNonZeros = matrix().nonZeros();

  auto weightsReal = std::vector<eckit::linalg::Triplet>(nNonZeros);
  auto weightsImag = std::vector<eckit::linalg::Triplet>(nNonZeros);

  const auto sourceLonLats = array::make_view<double, 2>(source_.lonlat());
  const auto targetLonLats = array::make_view<double, 2>(target_.lonlat());

  // Make complex weights (would be nice if we could have a complex matrix).
  spaceMatrixForEach(matrix(),
                     [&](auto&& weight, auto i, auto j, auto dataIdx) {
    const auto sourceLonLat =
        PointLonLat(sourceLonLats(j, 0), sourceLonLats(j, 1));
    const auto targetLonLat =
        PointLonLat(targetLonLats(i, 0), targetLonLats(i, 1));

    const auto alpha = util::greatCircleCourse(sourceLonLat, targetLonLat);

    auto deltaAlpha =
        (alpha.first - alpha.second) * util::Constants::degreesToRadians();

    weightsReal[dataIdx] = {i, j, weight * std::cos(deltaAlpha)};
    weightsImag[dataIdx] = {i, j, weight * std::sin(deltaAlpha)};
  });

  // Deal with slightly old fashioned Matrix interface
  const auto buildMatrix = [&](auto& matrix, const auto& weights) {
    auto tempMatrix = Matrix(nRows, nCols, weights);
    matrix.swap(tempMatrix);
  };

  buildMatrix(matrixReal_, weightsReal);
  buildMatrix(matrixImag_, weightsImag);
}

void ParallelTransport::print(std::ostream&) const { ATLAS_NOTIMPLEMENTED; }

void ParallelTransport::do_execute(const FieldSet& sourceFieldSet,
                                   FieldSet& targetFieldSet,
                                   Metadata& metadata) const {
  ATLAS_TRACE("atlas::interpolation::method::ParallelTransport::do_execute()");

  const auto nFields = sourceFieldSet.size();
  ATLAS_ASSERT(nFields == targetFieldSet.size());

  for (auto i = 0; i < sourceFieldSet.size(); ++i) {
    do_execute(sourceFieldSet[i], targetFieldSet[i], metadata);
  }
}

void ParallelTransport::do_execute(const Field& sourceField, Field& targetField,
                                   Metadata&) const {
  ATLAS_TRACE("atlas::interpolation::method::ParallelTransport::do_execute()");

  if (!(sourceField.variables() == 2 || sourceField.variables() == 3)) {

    auto metadata = Metadata();
    Method::do_execute(sourceField, targetField, metadata);

    return;
  }

  if (target_.size() == 0) {
    return;
  }

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
void ParallelTransport::interpolate_vector_field(const Field& sourceField,
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
void ParallelTransport::interpolate_vector_field(const Field& sourceField,
                                                 Field& targetField) const {

  const auto sourceView = array::make_view<Value, Rank>(sourceField);
  auto targetView = array::make_view<Value, Rank>(targetField);
  targetView.assign(0.);

  // Matrix multiplication split in two to simulate complex variable
  // multiplication.
  matrixMultiply(matrixReal_, sourceView, targetView,
                 [](const auto& weight, auto&& sourceVars, auto&& targetVars) {
    targetVars(0) += weight * sourceVars(0);
    targetVars(1) += weight * sourceVars(1);
  });
  matrixMultiply(matrixImag_, sourceView, targetView,
                 [](const auto& weight, auto&& sourceVars, auto&& targetVars) {
    targetVars(0) -= weight * sourceVars(1);
    targetVars(1) += weight * sourceVars(0);
  });

  if (sourceField.variables() == 3) {
    matrixMultiply(
        matrix(), sourceView, targetView,
        [](const auto& weight, auto&& sourceVars,
           auto&& targetVars) { targetVars(2) = weight * sourceVars(2); });
  }
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
