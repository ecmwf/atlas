/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/method/Method.h"

namespace atlas {
namespace interpolation {
namespace method {

/// @brief   Cubed sphere bilinear interpolation method.
///
/// @details Performs bilinear (quads) and baycentric (triangles) interpolation
///          accross cubed sphere tiles in (alpha, beta) coordinates.
///          Adding int "halo" (default 0) to the config controls the amount of
///          halo to consider when seraching for interpolation polygons. Adding
///          int "list_size" (defualt 8) will change the number of cells
///          returned by the internal kd-tree search. Increasing both values
///          may resolve errors if setup method fails to find cells.
///          The automatic halo exchange in the execute method can be disabled
///          by setting "halo_exchange" to false.
class CubedSphereBilinear : public Method {
public:
    CubedSphereBilinear(const Config& config): Method(config) {
        config.get("halo", halo_);
        config.get("list_size", listSize_);
        config.get("halo_exchange", halo_exchange_);
    }
    virtual ~CubedSphereBilinear() override {}

    virtual void print(std::ostream&) const override;
    virtual const FunctionSpace& source() const override { return source_; }
    virtual const FunctionSpace& target() const override { return target_; }

private:
    using Method::do_setup;
    void do_setup(const FunctionSpace& source, const FunctionSpace& target) override;
    void do_setup(const Grid& source, const Grid& target, const Cache&) override;
    void do_setup(const FunctionSpace& source, const FunctionSpace& target, const Cache&) override;

    FunctionSpace source_;
    FunctionSpace target_;

    int halo_{0};
    int listSize_{8};
    bool halo_exchange_{true};
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
