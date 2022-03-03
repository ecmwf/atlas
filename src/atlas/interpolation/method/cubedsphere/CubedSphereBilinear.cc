/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/interpolation/method/cubedsphere/CubedSphereBilinear.h"

namespace atlas {
namespace interpolation {
namespace method {

namespace  {
MethodBuilder<CubedSphereBilinear> __builder("cubedsphere-bilinear");
} // namespace

void CubedSphereBilinear::do_setup(const FunctionSpace &source, const FunctionSpace &target) {}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
