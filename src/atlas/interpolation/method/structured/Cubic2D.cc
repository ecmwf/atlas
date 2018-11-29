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
#include "kernels/CubicHorizontalKernel.h"


namespace atlas {
namespace interpolation {
namespace method {

using Cubic2D = StructuredInterpolation2D<CubicHorizontalKernel>;

namespace {

static MethodBuilder<Cubic2D> __builder1( "structured-cubic2D" );
static MethodBuilder<Cubic2D> __builder2( "cubic2D" );
static MethodBuilder<Cubic2D> __builder3( "structured-bicubic" );
static MethodBuilder<Cubic2D> __builder4( "bicubic" );

}  // namespace

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
