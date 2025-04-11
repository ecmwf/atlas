/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#pragma once

#include "StructuredInterpolation2D.h"
#include "kernels/CubicHorizontalKernel.h"

namespace atlas {
namespace interpolation {
namespace method {

class Cubic2D : public StructuredInterpolation2D<CubicHorizontalKernel> {
public:
    Cubic2D(const Config&);

    virtual std::string name() const override { return "structured-bicubic"; };
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
