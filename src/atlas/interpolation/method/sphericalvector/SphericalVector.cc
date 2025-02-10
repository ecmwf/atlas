/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <variant>

#include "atlas/interpolation/method/sphericalvector/SphericalVector.h"

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/interpolation/Cache.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/interpolation/method/sphericalvector/ComplexMatrixMultiply.h"
#include "atlas/interpolation/method/sphericalvector/Types.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/option/Options.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/PackVectorFields.h"
#include "eckit/config/LocalConfiguration.h"

namespace atlas {
namespace interpolation {
namespace method {

namespace {
MethodBuilder<SphericalVector> __builder("spherical-vector");
}

using namespace detail;

using WeightsMatMul = detail::ComplexMatrixMultiply<true>;
using WeightsMatMulAdjoint = detail::ComplexMatrixMultiply<false>;

SphericalVector::SphericalVector(const Config& config) : Method(config) {
  const auto* conf = dynamic_cast<const eckit::LocalConfiguration*>(&config);
  ATLAS_ASSERT(conf, "config must be derived from eckit::LocalConfiguration");
  interpolationScheme_ = conf->getSubConfiguration("scheme");
  adjoint_ = conf->getBool("adjoint", false);
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
  const auto m = atlas::linalg::make_non_owning_eckit_sparse_matrix(matrix());
  const auto nRows = static_cast<Index>(m.rows());
  const auto nCols = static_cast<Index>(m.cols());
  const auto nNonZeros = static_cast<std::size_t>(m.nonZeros());
  const auto* outerIndices = m.outer();
  const auto* innerIndices = m.inner();
  const auto* baseWeights  = m.data();

  // Note: need to store copy of weights as Eigen3 sorts compressed rows by j
  // whereas eckit does not.
  auto complexTriplets = ComplexTriplets(nNonZeros);
  auto realTriplets = RealTriplets(nNonZeros);

  const auto sourceLonLatsView = array::make_view<double, 2>(source_.lonlat());
  const auto targetLonLatsView = array::make_view<double, 2>(target_.lonlat());
  const auto unitSphere = geometry::UnitSphere{};

  atlas_omp_parallel_for(auto rowIndex = Index{0}; rowIndex < nRows;
                         ++rowIndex) {
    for (auto dataIndex = outerIndices[rowIndex];
         dataIndex < outerIndices[rowIndex + 1]; ++dataIndex) {
      const auto colIndex = static_cast<Index>(innerIndices[dataIndex]);
      const auto baseWeight = baseWeights[dataIndex];

      const auto sourceLonLat = PointLonLat(sourceLonLatsView(colIndex, 0),
                                            sourceLonLatsView(colIndex, 1));
      const auto targetLonLat = PointLonLat(targetLonLatsView(rowIndex, 0),
                                            targetLonLatsView(rowIndex, 1));

      const auto alpha =
          unitSphere.greatCircleCourse(sourceLonLat, targetLonLat);

      const auto deltaAlpha =
          (alpha.first - alpha.second) * util::Constants::degreesToRadians();

      complexTriplets[dataIndex] = ComplexTriplet{
          rowIndex, colIndex,
          baseWeight * Complex{std::cos(deltaAlpha), std::sin(deltaAlpha)}};
      realTriplets[dataIndex] = RealTriplet{rowIndex, colIndex, baseWeight};
    }
  }

  complexWeights_ = ComplexMatrix(nRows, nCols, complexTriplets);
  realWeights_ = RealMatrix(nRows, nCols, realTriplets);

  if (adjoint_) {
    complexWeightsAdjoint_ = complexWeights_.adjoint();
    realWeightsAdjoint_ = realWeights_.adjoint();
  }
}

void SphericalVector::print(std::ostream&) const { ATLAS_NOTIMPLEMENTED; }

void SphericalVector::do_execute(const FieldSet& sourceFieldSet,
                                 FieldSet& targetFieldSet,
                                 Metadata& metadata) const {
  ATLAS_TRACE("atlas::interpolation::method::SphericalVector::do_execute()");
  ATLAS_ASSERT(sourceFieldSet.size() == targetFieldSet.size());

  const auto packedSourceFieldSet = util::pack_vector_fields(sourceFieldSet);
  auto packedTargetFieldSet = util::pack_vector_fields(targetFieldSet);

  for (auto i = 0; i < packedSourceFieldSet.size(); ++i) {
    do_execute(packedSourceFieldSet[i], packedTargetFieldSet[i], metadata);
  }

  util::unpack_vector_fields(packedTargetFieldSet, targetFieldSet);
}

void SphericalVector::do_execute(const Field& sourceField, Field& targetField,
                                 Metadata&) const {
  ATLAS_TRACE("atlas::interpolation::method::SphericalVector::do_execute()");

  if (targetField.size() == 0) {
    return;
  }

  const auto fieldType = sourceField.metadata().getString("type", "");
  if (fieldType != "vector") {
    auto metadata = Metadata();
    Method::do_execute(sourceField, targetField, metadata);

    return;
  }

  Method::check_compatibility(sourceField, targetField, matrix());

  haloExchange(sourceField);
  interpolate_vector_field(sourceField, targetField,
                           WeightsMatMul(complexWeights_, realWeights_));
  targetField.set_dirty();
}

void SphericalVector::do_execute_adjoint(FieldSet& sourceFieldSet,
                                         const FieldSet& targetFieldSet,
                                         Metadata& metadata) const {
  ATLAS_TRACE(
      "atlas::interpolation::method::SphericalVector::do_execute_adjoint()");
  ATLAS_ASSERT(sourceFieldSet.size() == targetFieldSet.size());

  auto packedSourceFieldSet = util::pack_vector_fields(sourceFieldSet);
  const auto packedTargetFieldSet = util::pack_vector_fields(targetFieldSet);

  for (auto i = 0; i < packedSourceFieldSet.size(); ++i) {
    do_execute_adjoint(packedSourceFieldSet[i], packedTargetFieldSet[i],
                       metadata);
  }

  util::unpack_vector_fields(packedSourceFieldSet, sourceFieldSet);
}

void SphericalVector::do_execute_adjoint(Field& sourceField,
                                         const Field& targetField,
                                         Metadata& metadata) const {
  ATLAS_TRACE(
      "atlas::interpolation::method::SphericalVector::do_execute_adjoint()");

  if (targetField.size() == 0) {
    return;
  }

  const auto fieldType = sourceField.metadata().getString("type", "");
  if (fieldType != "vector") {
    auto metadata = Metadata();
    Method::do_execute_adjoint(sourceField, targetField, metadata);

    return;
  }

  Method::check_compatibility(sourceField, targetField, matrix());

  ATLAS_ASSERT(adjoint_, "\"adjoint\" needs to be set to \"true\" in Config.");
  interpolate_vector_field(
      targetField, sourceField,
      WeightsMatMulAdjoint(complexWeightsAdjoint_, realWeightsAdjoint_));
  adjointHaloExchange(sourceField);
}

template <typename MatMul>
void SphericalVector::interpolate_vector_field(const Field& sourceField,
                                               Field& targetField,
                                               const MatMul& matMul) {
  const auto sourceViewVariant = array::make_view_variant(sourceField);

  const auto sourceViewVisitor = [&](auto sourceView) {
    if constexpr (array::is_rank<2, 3>(sourceView) &&
                  array::is_non_const_value_type<float, double>(sourceView)) {
      using SourceView = std::decay_t<decltype(sourceView)>;
      using Value = typename SourceView::non_const_value_type;
      constexpr auto Rank = SourceView::rank();
      auto targetView = array::make_view<Value, Rank>(targetField);

      switch (sourceField.variables()) {
        case 2:
          return matMul.apply(sourceView, targetView, twoVector);
        case 3:
          return matMul.apply(sourceView, targetView, threeVector);
        default:
          ATLAS_THROW_EXCEPTION("Error: no support for " +
                                std::to_string(sourceField.variables()) +
                                " variable vector fields.\n" +
                                " Number of variables must be 2 or 3.");
      }

    } else {
      ATLAS_THROW_EXCEPTION(
          "Error: no support for rank = " + std::to_string(sourceField.rank()) +
          " and value type = " + sourceField.datatype().str() + ".\n" +
          "Vector field must have rank 2 or 3 with value type "
          "float or double");
    }
  };

  std::visit(sourceViewVisitor, sourceViewVariant);
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
