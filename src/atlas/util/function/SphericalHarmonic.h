/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include <functional>

namespace atlas {

namespace util {

namespace function {

/// \brief An analytic function that provides a spherical harmonic on a 2D sphere
///
/// \param n Total wave number
/// \param m Zonal wave number
/// \param lon Longitude in degrees
/// \param lat Latitude in degrees
/// \return spherical harmonic
double spherical_harmonic(int n, int m, double lon, double lat);

/// \brief SphericalHarmonic operator that provides a spherical harmonic on a 2D sphere
class SphericalHarmonic {
public:
    /// \brief Constructor
    ///
    /// \param n          Total wave number
    /// \param m          Zonal wave number
    /// \param caching    When true, internally cache the results of Associated Legendre Polynomials
    ///                   in a map. Warning: this is not thread-safe
    SphericalHarmonic(int n, int m, bool caching = false);

    /// \brief Evaluate the spherical harmonic function for given longitude and latitude
    double operator()(double lon, double lat) const;

private:
    std::function<double(double, double)> Y_;
};

}  // namespace function

}  // namespace util

}  // namespace atlas
