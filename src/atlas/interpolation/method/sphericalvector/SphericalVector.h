/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/library/defines.h"
#if ATLAS_HAVE_EIGEN

#pragma once

#include <complex>
#include <memory>

#include <Eigen/Sparse>

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/method/Method.h"
#include "atlas/linalg/sparse.h"
#include "eckit/config/Configuration.h"

namespace atlas {
namespace interpolation {
namespace method {

class SphericalVector : public Method {
 public:
  using Complex = std::complex<double>;

  template <typename Value>
  using SparseMatrix = Eigen::SparseMatrix<Value, Eigen::RowMajor>;
  using ComplexMatrix = SparseMatrix<Complex>;

  /// @brief   Interpolation post-processor for vector field interpolation
  ///
  /// @details Takes a base interpolation config keyed to "scheme" and creates
  ///          A set of complex intepolation weights which rotate source vector
  ///          elements into the target elements' individual frame of reference.
  ///          Method works by creating a great-circle arc between the source
  ///          and target locations; the amount of rotation is determined by the
  ///          difference the in the great-cricle course (relative to north) at
  ///          the source and target location.
  ///          Both source and target fields require the "type" metadata to be
  ///          set to "vector" for this method to be invoked. Otherwise, the
  ///          base scalar interpolation method is invoked.
  ///          Note: This method only works with matrix-based interpolation
  ///          schemes.
  ///
  SphericalVector(const Config& config) : Method(config) {
    const auto& conf = dynamic_cast<const eckit::LocalConfiguration&>(config);
    interpolationScheme_ = conf.getSubConfiguration("scheme");
  }
  virtual ~SphericalVector() override {}

  void print(std::ostream&) const override;
  const FunctionSpace& source() const override { return source_; }
  const FunctionSpace& target() const override { return target_; }

  void do_execute(const FieldSet& sourceFieldSet, FieldSet& targetFieldSet,
                  Metadata& metadata) const override;
  void do_execute(const Field& sourceField, Field& targetField,
                  Metadata& metadata) const override;

  void do_execute_adjoint(FieldSet& sourceFieldSet,
                          const FieldSet& targetFieldSet,
                          Metadata& metadata) const override {
    ATLAS_NOTIMPLEMENTED;
  }
  void do_execute_adjoint(Field& sourceField, const Field& targetField,
                          Metadata& metadata) const override {
    ATLAS_NOTIMPLEMENTED;
  }

 private:
  template <typename Value>
  void interpolate_vector_field(const Field& source, Field& target) const;

  template <typename Value, int Rank>
  void interpolate_vector_field(const Field& source, Field& target) const;

  using Method::do_setup;
  void do_setup(const FunctionSpace& source,
                const FunctionSpace& target) override;
  void do_setup(const Grid& source, const Grid& target, const Cache&) override;

  eckit::LocalConfiguration interpolationScheme_;

  FunctionSpace source_;
  FunctionSpace target_;

  std::shared_ptr<ComplexMatrix> complexWeights_;
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas

#endif