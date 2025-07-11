/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include "Cubic2D.h"

#include "atlas/interpolation/method/MethodFactory.h"

namespace atlas {
namespace interpolation {
namespace method {

namespace {

using util::FactoryDeprecated;

MethodBuilder<Cubic2D> __builder1("structured-bicubic");
MethodBuilder<Cubic2D> __builder2("structured-cubic2D", FactoryDeprecated("Please use structured-bicubic"));
MethodBuilder<Cubic2D> __builder3("cubic2D", FactoryDeprecated("Please use structured-bicubic"));
MethodBuilder<Cubic2D> __builder4("bicubic", FactoryDeprecated("Please use structured-bicubic"));

}  // namespace

Cubic2D::Cubic2D(const Config& config): StructuredInterpolation2D<CubicHorizontalKernel>(config) {}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
