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

#include "StructuredInterpolation3D.h"
#include "kernels/Cubic3DKernel.h"

namespace atlas {
namespace interpolation {
namespace method {

class Cubic3D : public StructuredInterpolation3D<Cubic3DKernel> {
public:
    Cubic3D(const Config&);
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
