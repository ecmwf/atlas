/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include "Linear3D.h"

#include "atlas/interpolation/method/MethodFactory.h"

namespace atlas {
namespace interpolation {
namespace method {

namespace {

MethodBuilder<Linear3D> __builder1("structured-linear3D");
MethodBuilder<Linear3D> __builder2("linear3D");
MethodBuilder<Linear3D> __builder3("structured-trilinear");
MethodBuilder<Linear3D> __builder4("trilinear");

}  // namespace

Linear3D::Linear3D(const Config& config): StructuredInterpolation3D<Linear3DKernel>(config) {}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
