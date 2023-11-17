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

  atlas_omp_parallel_for (auto i = 0; i < nRows; ++i) {
    for (auto dataIdx = rowIndices[i]; dataIdx < rowIndices[i+1]; ++dataIdx) {
      const auto j = colIndices[dataIdx];
      auto&& value = valData[dataIdx];
      functor(value, i, j);
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

  auto weights00 = std::vector<eckit::linalg::Triplet>();
  auto weights01 = std::vector<eckit::linalg::Triplet>();
  auto weights10 = std::vector<eckit::linalg::Triplet>();
  auto weights11 = std::vector<eckit::linalg::Triplet>();

  const auto sourceLonLats = array::make_view<double, 2>(source_.lonlat());
  const auto targetLonLats = array::make_view<double, 2>(target_.lonlat());


  spaceMatrixForEach(matrix(), [&](auto&& weight, int i, int j){
    const auto sourceLonLat = PointLonLat(sourceLonLats(j, 0), sourceLonLats(j, 1));
    const auto targetLonLat = PointLonLat(targetLonLats(i, 0), targetLonLats(i, 1));

    const auto alpha = util::greatCircleCourse(sourceLonLat, targetLonLat);

    auto deltaAlpha = (alpha.first - alpha.second) * util::Constants::degreesToRadians();


    weights00.emplace_back(i, j, weight * std::cos(deltaAlpha));
    weights01.emplace_back(i, j, -weight * std::sin(deltaAlpha));
    weights10.emplace_back(i, j, weight * std::sin(deltaAlpha));
    weights11.emplace_back(i, j, weight * std::cos(deltaAlpha));
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

void ParallelTransport::do_execute(const Field &source, Field &target, Metadata &) const
{
  ATLAS_TRACE("atlas::interpolation::method::ParallelTransport::do_execute()");

  if (!(source.variables() == 2 || source.variables() == 3)) {

    auto metadata = Metadata();
    Method::do_execute(source, target, metadata);

    return;

  }

  if (target_.size() == 0) {
    return;
  }

  haloExchange(source);

  if (source.datatype().kind() == array::DataType::KIND_REAL64) {
    interpolate_vector_field<double>(source, target);
  }
  else if (source.datatype().kind() == array::DataType::KIND_REAL32) {
    interpolate_vector_field<float>(source, target);
  }
  else {
    ATLAS_NOTIMPLEMENTED;
  }

  target.set_dirty();
}

template<typename Value>
void ParallelTransport::interpolate_vector_field(const Field &source, Field &target) const
{
  if (source.rank() == 2) {
     interpolate_vector_field<Value, 2>(source, target);
  }
  else if (source.rank() == 3 ) {
    interpolate_vector_field<Value, 3>(source, target);
  }
  else {
    ATLAS_NOTIMPLEMENTED;
  }
}

template<typename Value, int Rank>
void ParallelTransport::interpolate_vector_field(const Field &source, Field &target) const
{
  using namespace linalg;

  const auto sourceView = array::make_view<Value, Rank>(source);
  auto targetView = array::make_view<Value, Rank>(target);

  array::make_view<Value, Rank>(target).assign(0);



  const auto matrixMultiply = [&](const auto& matrix, int sourceVariableIdx, int targetVariableIdx) {

    spaceMatrixForEach(matrix, [&](const auto& weight, int i, int j){

      const auto adder = [&](auto&& targetElem, auto&& sourceElem){
        targetElem += weight * sourceElem;
      };

      if constexpr (Rank == 2) {
        adder(targetView(i, targetVariableIdx), sourceView(j, sourceVariableIdx));
      }
      else {
        const auto sourceSlice = sourceView.slice(j, array::Range::all(), sourceVariableIdx);
        auto targetSlice = targetView.slice(i, array::Range::all(), targetVariableIdx);
        array::helpers::ArrayForEach<0>::apply(std::tie(targetSlice, sourceSlice), adder);
      }

    });
  };

  matrixMultiply(matrix00_, 0, 0);
  matrixMultiply(matrix01_, 1, 0);
  matrixMultiply(matrix10_, 0, 1);
  matrixMultiply(matrix11_, 1, 1);

  if (source.variables() == 3) {
    matrixMultiply(matrix(), 2, 2);
  }

}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
