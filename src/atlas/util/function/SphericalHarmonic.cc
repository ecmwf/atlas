/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <map>

#include "atlas/runtime/Exception.h"
#include "atlas/util/Constants.h"
#include "atlas/util/function/SphericalHarmonic.h"

namespace atlas {

namespace util {

namespace function {


namespace {
static double factorial(double v) {
    if (v == 0) {
        return 1;
    }
    double result = v;
    while (--v > 0) {
        result *= v;
    }
    return result;
}

static double double_factorial(double x) {
    if (x == 0 || x == -1) {
        return 1;
    }

    double result = x;
    while ((x -= 2) > 0) {
        result *= x;
    }
    return result;
}

// Associated Legendre Polynomial
static double P(const int n, const int m, const double x) {
    // No recursive calculation needed
    if (n == m) {
        return (std::pow(-1.0, m) * double_factorial(2 * m - 1) * std::pow(std::sqrt(1. - x * x), m));
    }

    if (n == m + 1) {
        return x * (2 * m + 1) * P(m, m, x);
    }

    // Formula 1
    return (x * (2 * n - 1) * P(n - 1, m, x) - (n + m - 1) * P(n - 2, m, x)) / (n - m);
}

static double K(const int n, const int m) {
    //When m is less than 0, multiply - 1 to pass in
    return std::sqrt(((2 * n + 1) * factorial(n - m)) / (4 * M_PI * factorial(n + m)));
}
}  // namespace

double spherical_harmonic(int n, int m, double lon, double lat) {
    const int abs_m = std::abs(m);

    ATLAS_ASSERT(n >= abs_m);

    double colat = (90. - lat) * Constants::degreesToRadians();
    lon *= Constants::degreesToRadians();

    if (m == 0) {
        return (K(n, 0) * P(n, 0, std::cos(colat)));
    }

    if (m > 0) {
        return (M_SQRT2 * K(n, m) * std::cos(m * lon) * P(n, m, std::cos(colat)));
    }

    // When m is less than 0, multiply - 1 in advance and send it to K
    return (M_SQRT2 * K(n, abs_m) * std::sin(abs_m * lon) * P(n, abs_m, std::cos(colat)));
}

SphericalHarmonic::SphericalHarmonic(int n, int m, bool caching) {
    const int abs_m = std::abs(m);

    ATLAS_ASSERT(n >= abs_m);

    const double Knm                  = K(n, abs_m);
    std::function<double(double)> Pnm = [n, abs_m](double x) { return P(n, abs_m, x); };

    if (caching) {
        // Computation of associated legendre polynomials is
        // very expensive. With this trick (off by default),
        // we can cache repeated values, useful for structured
        // grids with constant latitude.
        Pnm = [Pnm](double x) {
            static std::map<double, double> memo;

            auto it = memo.find(x);
            if (it != memo.end()) {
                return it->second;
            }
            auto result = Pnm(x);
            memo[x]     = result;
            return result;
        };
    }

    if (m == 0) {
        Y_ = [Knm, Pnm](double lon, double colat) { return Knm * Pnm(std::cos(colat)); };
    }
    else if (m > 0) {
        Y_ = [m, Knm, Pnm](double lon, double colat) {
            return M_SQRT2 * Knm * std::cos(m * lon) * Pnm(std::cos(colat));
        };
    }
    else {
        Y_ = [m, Knm, Pnm](double lon, double colat) {
            return M_SQRT2 * Knm * std::sin(-m * lon) * Pnm(std::cos(colat));
        };
    }
}

double SphericalHarmonic::operator()(double lon, double lat) const {
    double colat = (90. - lat) * Constants::degreesToRadians();
    lon *= Constants::degreesToRadians();
    return Y_(lon, colat);
}


}  // namespace function

}  // namespace util

}  // namespace atlas
