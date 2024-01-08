/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/sphericalvector/SphericalVector.h"

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/interpolation/Cache.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/interpolation/method/sphericalvector/Types.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Geometry.h"
#include "eckit/config/LocalConfiguration.h"

namespace atlas {
namespace interpolation {
namespace method {

namespace  {
MethodBuilder<SphericalVector> __builder("spherical-vector");
}

using ComplexTriplets = detail::ComplexMatrix::Triplets;
using RealTriplets = detail::RealMatrix::Triplets;

SphericalVector::SphericalVector(const Config& config) : Method(config) {
  const auto& conf = dynamic_cast<const eckit::LocalConfiguration&>(config);
  interpolationScheme_ = conf.getSubConfiguration("scheme");
  adjoint_ = conf.getBool("adjoint", false);
}

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

  setMatrix(Interpolation(interpolationScheme_, source_, target_));

  // Get matrix data.
  const auto nRows = matrix().rows();
  const auto nCols = matrix().cols();
  const auto nNonZeros = matrix().nonZeros();
  const auto* outerIndices = matrix().outer();
  const auto* innerIndices = matrix().inner();
  const auto* baseWeights = matrix().data();

  using Index = detail::Index;

  // Note: need to store copy of weights as Eigen3 sorts compressed rows by j
  // whereas eckit does not.
  auto complexTriplets = ComplexTriplets(nNonZeros);
  auto realTriplets = RealTriplets(nNonZeros);

  // Make sure halo lonlats are same as owned points.
  auto sourceLonLats = source_.createField<double>(option::name("lonlat") |
                                                   option::variables(2));
  auto targetLonLats = target_.createField<double>(option::name("lonlat") |
                                                   option::variables(2));
  sourceLonLats.array().copy(source_.lonlat());
  targetLonLats.array().copy(target_.lonlat());
  sourceLonLats.haloExchange();
  targetLonLats.haloExchange();

  const auto sourceLonLatsView = array::make_view<double, 2>(sourceLonLats);
  const auto targetLonLatsView = array::make_view<double, 2>(targetLonLats);

  const auto unitSphere = geometry::UnitSphere{};

  atlas_omp_parallel_for(auto rowIndex = Index{0}; rowIndex < nRows;
                         ++rowIndex) {
    for (auto dataIndex = outerIndices[rowIndex];
         dataIndex < outerIndices[rowIndex + 1]; ++dataIndex) {
      const auto colIndex = innerIndices[dataIndex];
      const auto baseWeight = baseWeights[dataIndex];

      const auto sourceLonLat = PointLonLat(sourceLonLatsView(colIndex, 0),
                                            sourceLonLatsView(colIndex, 1));
      const auto targetLonLat = PointLonLat(targetLonLatsView(rowIndex, 0),
                                            targetLonLatsView(rowIndex, 1));

      const auto alpha =
          unitSphere.greatCircleCourse(sourceLonLat, targetLonLat);

      const auto deltaAlpha =
          (alpha.first - alpha.second) * util::Constants::degreesToRadians();

      complexTriplets[dataIndex] = {rowIndex, colIndex,
                                    {baseWeight * std::cos(deltaAlpha),
                                     baseWeight * std::sin(deltaAlpha)}};
      realTriplets[dataIndex] = {rowIndex, colIndex, baseWeight};
    }
  }

  const auto complexWeights =
      std::make_shared<detail::ComplexMatrix>(nRows, nCols, complexTriplets);

  const auto realWeights =
      std::make_shared<detail::RealMatrix>(nRows, nCols, realTriplets);

  weightsMatMul_= WeightsMatMul(complexWeights, realWeights);

  if (adjoint_) {

    const auto complexWeightsAdjoint =
        std::make_shared<detail::ComplexMatrix>(complexWeights->adjoint());

    const auto realWeightsAdjoint =
        std::make_shared<detail::RealMatrix>(realWeights->adjoint());

    weightsMatMulAdjoint_ =
        WeightsMatMulAdjoint(complexWeightsAdjoint, realWeightsAdjoint);
  }
}

void SphericalVector::print(std::ostream&) const { ATLAS_NOTIMPLEMENTED; }

void SphericalVector::do_execute(const FieldSet& sourceFieldSet,
                                 FieldSet& targetFieldSet,
                                 Metadata& metadata) const {
  ATLAS_TRACE("atlas::interpolation::method::SphericalVector::do_execute()");
  ATLAS_ASSERT(sourceFieldSet.size() == targetFieldSet.size());

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

  Method::check_compatibility(sourceField, targetField, matrix());

  haloExchange(sourceField);
  interpolate_vector_field(sourceField, targetField, weightsMatMul_);
  targetField.set_dirty();
}

void SphericalVector::do_execute_adjoint(FieldSet& sourceFieldSet,
                                         const FieldSet& targetFieldSet,
                                         Metadata& metadata) const {
  ATLAS_TRACE(
      "atlas::interpolation::method::SphericalVector::do_execute_adjoint()");
  ATLAS_ASSERT(sourceFieldSet.size() == targetFieldSet.size());

  for (auto i = 0; i < sourceFieldSet.size(); ++i) {
    do_execute_adjoint(sourceFieldSet[i], targetFieldSet[i], metadata);
  }
}

void SphericalVector::do_execute_adjoint(Field& sourceField,
                                         const Field& targetField,
                                         Metadata& metadata) const {
  ATLAS_TRACE(
      "atlas::interpolation::method::SphericalVector::do_execute_adjoint()");

  const auto fieldType = sourceField.metadata().getString("type", "");
  if (fieldType != "vector") {

    auto metadata = Metadata();
    Method::do_execute_adjoint(sourceField, targetField, metadata);

    return;
  }

  Method::check_compatibility(sourceField, targetField, matrix());

  ATLAS_ASSERT(adjoint_, "\"adjoint\" needs to be set to \"true\" in Config.");
  interpolate_vector_field(targetField, sourceField, weightsMatMulAdjoint_);
    adjointHaloExchange(sourceField);
}

template <typename MatMul>
void SphericalVector::interpolate_vector_field(const Field& sourceField,
                                               Field& targetField,
                                               const MatMul& matMul) {
  if (targetField.size() == 0) {
    return;
  }

  ATLAS_ASSERT_MSG(sourceField.variables() == 2 || sourceField.variables() == 3,
                   "Vector field can only have 2 or 3 components.");

  if (sourceField.datatype().kind() == array::DataType::KIND_REAL64) {
    interpolate_vector_field<double>(sourceField, targetField, matMul);
    return;
  }

  if (sourceField.datatype().kind() == array::DataType::KIND_REAL32) {
    interpolate_vector_field<float>(sourceField, targetField, matMul);
    return;
  }

  ATLAS_NOTIMPLEMENTED;
};

template <typename Value, typename MatMul>
void SphericalVector::interpolate_vector_field(const Field& sourceField,
                                               Field& targetField,
                                               const MatMul& matMul) {
  if (sourceField.rank() == 2) {
    interpolate_vector_field<Value, 2>(sourceField, targetField, matMul);
    return;
  }

  if (sourceField.rank() == 3) {
    interpolate_vector_field<Value, 3>(sourceField, targetField, matMul);
    return;
  }

  ATLAS_NOTIMPLEMENTED;
}

template <typename Value, int Rank, typename MatMul>
void SphericalVector::interpolate_vector_field(const Field& sourceField,
                                               Field& targetField,
                                               const MatMul& matMul) {
  const auto sourceView = array::make_view<Value, Rank>(sourceField);
  auto targetView = array::make_view<Value, Rank>(targetField);

  if (sourceField.variables() == 2) {
    matMul.apply(sourceView, targetView, detail::twoVector);
    return;
  }

  if (sourceField.variables() == 3) {
    matMul.apply(sourceView, targetView, detail::threeVector);
    return;
  }

  ATLAS_NOTIMPLEMENTED;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
