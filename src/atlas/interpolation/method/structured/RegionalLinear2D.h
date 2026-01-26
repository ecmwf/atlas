/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/interpolation/method/Method.h"

#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid/Grid.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

namespace atlas {
namespace interpolation {
namespace method {


class RegionalLinear2D : public Method {
 public:
  /// @brief   Regional linear interpolation
  ///
  /// @details 
  ///
  RegionalLinear2D(const Config& config) : Method(config), comm_(eckit::mpi::comm()) {}
  ~RegionalLinear2D() override {}

  void print(std::ostream&) const override;
  const FunctionSpace& source() const override { return source_; }
  const FunctionSpace& target() const override { return target_; }

  void do_execute(const FieldSet& sourceFieldSet, FieldSet& targetFieldSet,
                  Metadata& metadata) const override;
  void do_execute(const Field& sourceField, Field& targetField,
                  Metadata& metadata) const override;

  void do_execute_adjoint(FieldSet& sourceFieldSet,
                          const FieldSet& targetFieldSet,
                          Metadata& metadata) const override;
  void do_execute_adjoint(Field& sourceField, const Field& targetField,
                          Metadata& metadata) const override;

 private:
  using Method::do_setup;
  void do_setup(const FunctionSpace& source, const FunctionSpace& target) override;
  void do_setup(const Grid& source, const Grid& target, const Cache&) override;
  void do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) override;

  FunctionSpace source_{};
  FunctionSpace target_{};

  const eckit::mpi::Comm & comm_;
  size_t sourceSize_;
  std::vector<int> mpiTask_;
  size_t targetSize_;
  size_t sourceSendSize_;
  size_t targetRecvSize_;
  std::vector<int> sourceSendCounts_;
  std::vector<int> sourceSendDispls_;
  std::vector<int> targetRecvCounts_;
  std::vector<int> targetRecvDispls_;
  std::vector<size_t> sourceSendMapping_;
  std::vector<std::array<size_t,4>> stencil_;
  std::vector<std::array<double,4>> weights_;
  std::vector<size_t> stencilSize_;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
