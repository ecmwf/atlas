/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include "StructuredInterpolation2D.h"
#include "kernels/QuasiCubicHorizontalKernel.h"


namespace atlas {
namespace interpolation {
namespace method {

using QuasiCubic2D = StructuredInterpolation2D<QuasiCubicHorizontalKernel>;

namespace {

static MethodBuilder<QuasiCubic2D> __builder1( "structured-quasicubic2D" );
static MethodBuilder<QuasiCubic2D> __builder2( "quasicubic2D" );
static MethodBuilder<QuasiCubic2D> __builder3( "structured-biquasicubic" );
static MethodBuilder<QuasiCubic2D> __builder4( "biquasicubic" );

}  // namespace

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
