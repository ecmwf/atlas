/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include "StructuredInterpolation3D.h"
#include "kernels/Linear3DKernel.h"


namespace atlas {
namespace interpolation {
namespace method {

using Linear3D = StructuredInterpolation3D<Linear3DKernel>;

namespace {

static MethodBuilder<Linear3D> __builder1( "structured-linear3D" );
static MethodBuilder<Linear3D> __builder2( "linear3D" );
static MethodBuilder<Linear3D> __builder3( "structured-trilinear" );
static MethodBuilder<Linear3D> __builder4( "trilinear" );

}  // namespace

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
