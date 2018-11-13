/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include <cmath>
#include <limits>

#include "StructuredInterpolation.h"
#include "kernels/BicubicKernel.h"


namespace atlas {
namespace interpolation {
namespace method {

using Bicubic = StructuredInterpolation<BicubicKernel>;

namespace {

static MethodBuilder<Bicubic> __builder1( "structured-bicubic" );
static MethodBuilder<Bicubic> __builder2( "bicubic" );

}  // namespace

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
