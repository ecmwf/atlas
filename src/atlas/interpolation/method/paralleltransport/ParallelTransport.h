/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/method/Method.h"
#include "atlas/linalg/sparse.h"
#include "eckit/config/Configuration.h"

namespace atlas {
namespace interpolation {
namespace method {

class ParallelTransport : public Method {
 public:
  ParallelTransport(const Config& config): Method(config) {
    const auto& conf = dynamic_cast<const eckit::LocalConfiguration&>(config);
    interpolationScheme_ = conf.getSubConfiguration("scheme");

  }
  virtual ~ParallelTransport() override {}

  void print(std::ostream&) const override;
  const FunctionSpace& source() const override { return source_; }
  const FunctionSpace& target() const override { return target_; }

  void do_execute(const FieldSet& sourceFieldSet, FieldSet& targetFieldSet, Metadata& metadata) const override;
  void do_execute(const Field& sourceField, Field& targetField, Metadata& metadata) const override;

 private:

  template<typename Value>
  void interpolate_vector_field(const Field& source, Field& target) const;

  template<typename Value, int Rank>
  void interpolate_vector_field(const Field& source, Field& target) const;


  using Method::do_setup;
  void do_setup(const FunctionSpace& source, const FunctionSpace& target) override;
  void do_setup(const Grid& source, const Grid& target, const Cache&) override;

  eckit::LocalConfiguration interpolationScheme_;

  FunctionSpace source_;
  FunctionSpace target_;

  // Complex interpolation weights. We treat a (u, v) pair as a complex variable
  // with u on the real number line and v on the imaginary number line.
  Matrix matrixReal_;
  Matrix matrixImag_;

};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
