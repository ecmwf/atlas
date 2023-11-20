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

template<typename MatrixT, typename Functor>
void spaceMatrixForEach(MatrixT&& matrix, const Functor& functor) {

  const auto nRows = matrix.rows();
  const auto nCols = matrix.cols();
  const auto rowIndices = matrix.outer();
  const auto colIndices = matrix.inner();
  auto valData = matrix.data();

  atlas_omp_parallel_for (auto i = size_t{}; i < nRows; ++i) {
    for (auto dataIdx = rowIndices[i]; dataIdx < rowIndices[i+1]; ++dataIdx) {
      const auto j = size_t(colIndices[dataIdx]);
      auto& value = valData[dataIdx];

      if constexpr (std::is_invocable_v<Functor, decltype(value), size_t, size_t>) {
        functor(value, i, j);
      }
      else if constexpr (std::is_invocable_v<Functor, decltype(value), size_t, size_t, size_t>) {
        functor(value, i, j, dataIdx);
      }
      else {
        ATLAS_NOTIMPLEMENTED;
      }
    }
  }
}

}  // namespace

void ParallelTransport::do_setup(const Grid& source, const Grid& target, const Cache&) {
  ATLAS_NOTIMPLEMENTED;
}

void ParallelTransport::do_setup(const FunctionSpace& source, const FunctionSpace& target) {
  ATLAS_TRACE("interpolation::method::ParallelTransport::do_setup");
  source_ = source;
  target_ = target;

  if (target_.size() == 0) {
    return;
  }

  const auto baseInterpolator = Interpolation(interpolationScheme_, source_, target_);
  setMatrix(MatrixCache(baseInterpolator));

  // Get matrix dimensions.
  const auto nRows = matrix().rows();
  const auto nCols = matrix().cols();
  const auto nNonZeros = matrix().nonZeros();

  auto weights00 = std::vector<eckit::linalg::Triplet>(nNonZeros);
  auto weights01 = std::vector<eckit::linalg::Triplet>(nNonZeros);
  auto weights10 = std::vector<eckit::linalg::Triplet>(nNonZeros);
  auto weights11 = std::vector<eckit::linalg::Triplet>(nNonZeros);

  const auto sourceLonLats = array::make_view<double, 2>(source_.lonlat());
  const auto targetLonLats = array::make_view<double, 2>(target_.lonlat());


  spaceMatrixForEach(matrix(), [&](auto&& weight, auto i, auto j, auto dataIdx){
    const auto sourceLonLat = PointLonLat(sourceLonLats(j, 0), sourceLonLats(j, 1));
    const auto targetLonLat = PointLonLat(targetLonLats(i, 0), targetLonLats(i, 1));

    const auto alpha = util::greatCircleCourse(sourceLonLat, targetLonLat);

    auto deltaAlpha = (alpha.first - alpha.second) * util::Constants::degreesToRadians();

    weights00[dataIdx] = {i, j, weight * std::cos(deltaAlpha)};
    weights01[dataIdx] = {i, j, -weight * std::sin(deltaAlpha)};
    weights10[dataIdx] = {i, j, weight * std::sin(deltaAlpha)};
    weights11[dataIdx] = {i, j, weight * std::cos(deltaAlpha)};
  });


  // Deal with slightly old fashioned Matrix interface
  const auto buildMatrix = [&](auto& matrix, const auto& weights){
    auto tempMatrix = Matrix(nRows, nCols, weights);
    matrix.swap(tempMatrix);
  };

  buildMatrix(matrix00_, weights00);
  buildMatrix(matrix10_, weights10);
  buildMatrix(matrix01_, weights01);
  buildMatrix(matrix11_, weights11);

}

void ParallelTransport::print(std::ostream&) const {
  ATLAS_NOTIMPLEMENTED;
}

void ParallelTransport::do_execute(const FieldSet& sourceFieldSet, FieldSet& targetFieldSet, Metadata& metadata) const
{
  ATLAS_TRACE("atlas::interpolation::method::ParallelTransport::do_execute()");

  const auto nFields = sourceFieldSet.size();
  ATLAS_ASSERT(nFields == targetFieldSet.size());

  for (auto i = 0; i < sourceFieldSet.size(); ++i) {
    do_execute(sourceFieldSet[i], targetFieldSet[i], metadata);
  }
}

void ParallelTransport::do_execute(const Field& sourceField, Field& targetField, Metadata &) const
{
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
  }
  else if (sourceField.datatype().kind() == array::DataType::KIND_REAL32) {
    interpolate_vector_field<float>(sourceField, targetField);
  }
  else {
    ATLAS_NOTIMPLEMENTED;
  }

  targetField.set_dirty();
}

template<typename Value>
void ParallelTransport::interpolate_vector_field(const Field& sourceField, Field& targetField) const
{
  if (sourceField.rank() == 2) {
     interpolate_vector_field<Value, 2>(sourceField, targetField);
  }
  else if (sourceField.rank() == 3 ) {
    interpolate_vector_field<Value, 3>(sourceField, targetField);
  }
  else {
    ATLAS_NOTIMPLEMENTED;
  }
}

template<typename Value, int Rank>
void ParallelTransport::interpolate_vector_field(const Field& sourceField, Field& targetField) const
{
  using namespace linalg;

  const auto sourceView = array::make_view<Value, Rank>(sourceField);
  auto targetView = array::make_view<Value, Rank>(targetField);

  array::make_view<Value, Rank>(targetField).assign(0);



  const auto matrixMultiply = [&](const auto& matrix, auto sourceVariableIdx, auto targetVariableIdx) {

    spaceMatrixForEach(matrix, [&](const auto& weight, auto i, auto j){

      const auto adder = [&](auto&& targetElem, auto&& sourceElem){
        targetElem += weight * sourceElem;
      };

      if constexpr (Rank == 2) {
        adder(targetView(i, targetVariableIdx), sourceView(j, sourceVariableIdx));
      }
      else if constexpr (Rank == 3) {
        const auto sourceSlice = sourceView.slice(j, array::Range::all(), sourceVariableIdx);
        auto targetSlice = targetView.slice(i, array::Range::all(), targetVariableIdx);
        array::helpers::ArrayForEach<0>::apply(std::tie(targetSlice, sourceSlice), adder);
      }
      else {
        ATLAS_NOTIMPLEMENTED;
      }

    });
  };

  matrixMultiply(matrix00_, 0, 0);
  matrixMultiply(matrix01_, 1, 0);
  matrixMultiply(matrix10_, 0, 1);
  matrixMultiply(matrix11_, 1, 1);

  if (sourceField.variables() == 3) {
    matrixMultiply(matrix(), 2, 2);
  }

}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
