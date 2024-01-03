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

  using Index = std::decay_t<decltype(*innerIndices)>;

  // Note: need to store copy of weights as Eigen3 sorts compressed rows by j
  // whereas eckit does not.
  auto complexTriplets = ComplexTriplets(nNonZeros);
  auto realTriplets = RealTriplets(nNonZeros);

  const auto sourceLonLats = array::make_view<double, 2>(source_.lonlat());
  const auto targetLonLats = array::make_view<double, 2>(target_.lonlat());

  const auto unitSphere = geometry::UnitSphere{};

  atlas_omp_parallel_for(auto rowIndex = Index{0}; rowIndex < nRows;
                         ++rowIndex) {
    for (auto dataIndex = outerIndices[rowIndex];
         dataIndex < outerIndices[rowIndex + 1]; ++dataIndex) {
      const auto colIndex = innerIndices[dataIndex];
      const auto baseWeight = baseWeights[dataIndex];

      const auto sourceLonLat =
          PointLonLat(sourceLonLats(colIndex, 0), sourceLonLats(colIndex, 1));
      const auto targetLonLat =
          PointLonLat(targetLonLats(rowIndex, 0), targetLonLats(rowIndex, 1));

      const auto alpha =
          unitSphere.greatCircleCourse(sourceLonLat, targetLonLat);

      const auto deltaAlpha =
          (alpha.first - alpha.second) * util::Constants::degreesToRadians();

      complexTriplets[dataIndex] = {rowIndex, colIndex,
                                    std::polar(baseWeight, deltaAlpha)};
      realTriplets[dataIndex] = {rowIndex, colIndex, baseWeight};
    }
  }

  const auto complexWeights =
      std::make_shared<detail::ComplexMatrix>(nRows, nCols, complexTriplets);

  const auto realWeights =
      std::make_shared<detail::RealMatrix>(nRows, nCols, realTriplets);

  weightsMatMul_= WeightsMatMul(complexWeights, realWeights);

}

SphericalVector::SphericalVector(const Config& config) : Method(config) {
  const auto& conf = dynamic_cast<const eckit::LocalConfiguration&>(config);
  interpolationScheme_ = conf.getSubConfiguration("scheme");
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

void SphericalVector::do_execute_adjoint(FieldSet& sourceFieldSet,
                                         const FieldSet& targetFieldSet,
                                         Metadata& metadata) const {
  ATLAS_NOTIMPLEMENTED;
}

void SphericalVector::do_execute_adjoint(Field& sourceField,
                                         const Field& targetField,
                                         Metadata& metadata) const {
  ATLAS_NOTIMPLEMENTED;
}

template <typename Value>
void SphericalVector::interpolate_vector_field(const Field& sourceField,
                                               Field& targetField) const {
  if (sourceField.rank() == 2) {
    interpolate_vector_field<Value, 2>(sourceField, targetField);
    return;
  }

  if (sourceField.rank() == 3) {
    interpolate_vector_field<Value, 3>(sourceField, targetField);
    return;
  }

  ATLAS_NOTIMPLEMENTED;
}

template <typename Value, int Rank>
void SphericalVector::interpolate_vector_field(const Field& sourceField,
                                               Field& targetField) const {
  const auto sourceView = array::make_view<Value, Rank>(sourceField);
  auto targetView = array::make_view<Value, Rank>(targetField);

  if (sourceField.variables() == 2) {
    weightsMatMul_.apply(sourceView, targetView, detail::twoVector);
    return;
  }

  if (sourceField.variables() == 3) {
    weightsMatMul_.apply(sourceView, targetView, detail::threeVector);
    return;
  }

  ATLAS_NOTIMPLEMENTED;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
