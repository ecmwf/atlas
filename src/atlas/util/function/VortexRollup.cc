/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/Constants.h"
#include "atlas/util/Earth.h"

#include "atlas/util/function/VortexRollup.h"

namespace atlas {

namespace util {

namespace function {

double vortex_rollup(double lon, double lat, double t) {
    lon *= Constants::degreesToRadians();
    lat *= Constants::degreesToRadians();

    auto sqr                  = [](const double x) { return x * x; };
    auto sech                 = [](const double x) { return 1. / std::cosh(x); };
    constexpr double two_pi   = 2. * M_PI;
    const double lambda_prime = std::atan2(-std::cos(lon - two_pi * t), std::tan(lat));
    const double rho          = 3. * std::sqrt(1. - sqr(std::cos(lat)) * sqr(std::sin(lon - two_pi * t)));
    double omega              = 0.;
    double a                  = Earth::radius();
    if (rho != 0.) {
        omega = 0.5 * 3 * std::sqrt(3) * a * two_pi * sqr(sech(rho)) * std::tanh(rho) / rho;
    }
    return -std::tanh(0.2 * rho * std::sin(lambda_prime - omega / a * t));
}

}  // namespace function

}  // namespace util

}  // namespace atlas
