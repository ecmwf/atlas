/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

namespace atlas {

namespace util {

namespace function {

/// \brief An analytic function that provides solid body rotation winds on a sphere.
///
/// All angles must be provided in degrees.
///
class SolidBodyRotation {
public:

    SolidBodyRotation() : SolidBodyRotation(0., 1.) {}


    SolidBodyRotation(const double beta) : SolidBodyRotation(beta, 1.) {}
    SolidBodyRotation(const double beta, const double radius);

    void wind(const double lon, const double lat, double& u, double& v) const;
    void vordiv(const double lon, const double lat, double& vor, double& div) const;
    double windMagnitude(const double lon, const double lat) const;
    double u(const double lon, const double lat) const;
    double v(const double lon, const double lat) const;
    double vorticity(const double lon, const double lat) const;
    double divergence(const double lon, const double lat) const;

    double windMagnitudeSquared(const double lon, const double lat) const;
    void windMagnitudeSquaredGradient(const double lon, const double lat, double& dfdx, double& dfdy) const;

private:
    double sin_beta_;
    double cos_beta_;
    double radius_;
};

}  // namespace function

}  // namespace util

}  // namespace atlas
