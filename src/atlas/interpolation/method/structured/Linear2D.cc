/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include "Linear2D.h"

#include "atlas/interpolation/method/MethodFactory.h"

namespace atlas {
namespace interpolation {
namespace method {

namespace {

using util::FactoryDeprecated;

MethodBuilder<Linear2D> __builder0("structured-linear");
MethodBuilder<Linear2D> __builder1("structured-bilinear", FactoryDeprecated("Please use structured-linear"));
MethodBuilder<Linear2D> __builder2("structured-linear2D", FactoryDeprecated("Please use structured-linear"));
MethodBuilder<Linear2D> __builder3("linear2D", FactoryDeprecated("Please use structured-linear"));
MethodBuilder<Linear2D> __builder4("bilinear", FactoryDeprecated("Please use structured-linear"));

}  // namespace

Linear2D::Linear2D(const Config& config): StructuredInterpolation2D<LinearHorizontalKernel>(config) {}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
