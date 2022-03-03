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

class CubedSphereBilinear : public Method {
public:
    CubedSphereBilinear(const Config& config) : Method(config) {}
    virtual ~CubedSphereBilinear() override {}

protected:
    virtual const FunctionSpace& source() const override { return source_; }
    virtual const FunctionSpace& target() const override { return target_; }

private:

    /**
   * @brief Create an interpolant sparse matrix relating two (pre-partitioned)
   * meshes,
   * using nearest neighbour method
   * @param source functionspace containing source elements
   * @param target functionspace containing target points
   */
    virtual void do_setup(const FunctionSpace& source, const FunctionSpace& target) override;
    virtual void do_setup(const Grid& source, const Grid& target, const Cache&) override;

    FunctionSpace source_;
    FunctionSpace target_;


};

}
} // namespace interpolation
} // namespace atlas

