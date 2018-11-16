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
#include "kernels/TricubicKernel.h"


namespace atlas {
namespace interpolation {
namespace method {

using Tricubic = StructuredInterpolation3D<TricubicKernel>;

namespace {

static MethodBuilder<Tricubic> __builder1( "structured-tricubic" );
static MethodBuilder<Tricubic> __builder2( "tricubic" );

}  // namespace

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
