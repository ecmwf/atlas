/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/Cache.h"
#include "atlas/interpolation/method/Method.h"
#include "atlas/grid/Grid.h"

#include "eckit/config/Configuration.h"


namespace atlas {
namespace interpolation {
namespace method {


class Binning : public Method {
 public:
  /// @brief   Binning procedure to carry out a regridding from high
  ///          to low resolution
  ///
  /// @details This method is part of the family of the Interpolation operations;
  ///          it relies on the evaluation of the transpose of an interpolation matrix.
  ///          The rigridding from high to low resolution is carried out by using
  ///          a 'binning matrix':
  ///                         binning matrix: B = N A^T W
  ///                    area weights matrix: W
  ///                   interpolation matrix: A
  ///            normalization factor matrix: N
  ///
  ///          Setup, configuration variables:
  ///                      <type>: method used to evaluate the 'B' matrix;
  ///                       value: 'binning'
  ///                    <scheme>: method used to evaluate the 'A' matrix;
  ///                       value: 'cubedsphere-bilinear', 'structured-bilinear', ...
  ///             <halo_exchange>: flag to control the halo exchange procedure
  ///                       value: 'true', 'false'
  ///                   <adjoint>: flag to control the adjoint operation
  ///                       value: 'true', 'false'
  ///
  Binning(const Config& config);
  ~Binning() override {}

  void print(std::ostream&) const override;
  const FunctionSpace& source() const override { return source_; }
  const FunctionSpace& target() const override { return target_; }

 private:
  using Method::do_setup;
  void do_setup(const FunctionSpace& source,
                const FunctionSpace& target) override;
  void do_setup(const Grid& source, const Grid& target, const Cache&) override;
  void do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) override;

  std::vector<double> getAreaWeights(const FunctionSpace& source) const;

  eckit::LocalConfiguration interpAncillaryScheme_{};

  FunctionSpace source_{};
  FunctionSpace target_{};
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
