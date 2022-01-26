/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/VortexRollup.h"

namespace atlas {

namespace util {

double vortex_rollup(double lon, double lat, double t,
                     const double mean, const bool degreesInput,
                     const double timescale) {

  if (degreesInput) {
    lon *= Constants::degreesToRadians();
    lat *= Constants::degreesToRadians();
  }

  auto sqr           = [](const double x) { return x * x; };
  auto sech          = [](const double x) { return 1. / std::cosh(x); };
  const double Omega = 2. * M_PI / timescale;
  t *= timescale;
  const double lambda_prime = std::atan2(-std::cos(lon - Omega * t), std::tan(lat));
  const double rho = 3. * std::sqrt(1. - sqr(std::cos(lat)) * sqr(std::sin(lon - Omega * t)));
  double omega     = 0.;
  double a         = Earth::radius();
  if (rho != 0.) {
    omega = 0.5 * 3 * std::sqrt(3) * a * Omega * sqr(sech(rho)) * std::tanh(rho) / rho;
  }
  double q = mean - std::tanh(0.2 * rho * std::sin(lambda_prime - omega / a * t));
  return q;
}

}

}
