/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>

#include "atlas/util/Constants.h"
#include "atlas/util/Earth.h"

#include "atlas/util/function/MDPI_functions.h"

namespace atlas {

namespace util {

namespace function {

double MDPI_sinusoid(double lon, double lat) {
    lon *= Constants::degreesToRadians();
    lat *= Constants::degreesToRadians();
    constexpr double length = 1.2 * M_PI;

    return 2. - std::cos(M_PI * std::acos(std::cos(lon) * std::cos(lat)) / length);
}

double MDPI_harmonic(double lon, double lat) {
    lon *= Constants::degreesToRadians();
    lat *= Constants::degreesToRadians();
    
    return 2. + std::pow(std::sin(2. * lat), 16) * std::cos(16. * lon);
}

double MDPI_vortex(double lon, double lat) {
    lon *= Constants::degreesToRadians();
    lat *= Constants::degreesToRadians();
    constexpr double dLon0  = 5.5;
    constexpr double dLat0  = 0.2;
    constexpr double dR0    = 3.0;
    constexpr double dD     = 5.0;
    constexpr double dT     = 6.0;
    constexpr double length = 1.2 * M_PI;

    auto sqr                  = [](const double x) { return x * x; };
    auto sech                 = [](const double x) { return 1. / std::cosh(x); };

    double dSinC  = std::sin(dLat0);
    double dCosC  = std::cos(dLat0);
    double dCosT = std::cos(lat);
    double dSinT = std::sin(lat);
    double dTrm = dCosT * std::cos(lon - dLon0);
    double dX = dSinC * dTrm - dCosC * dSinT;
    double dY = dCosT * std::sin(lon - dLon0);
    double dZ = dSinC * dSinT + dCosC * dTrm;

    double dlon = std::atan2(dY, dX);
    double dlat = std::asin(dZ);

    double dRho = dR0 * std::cos(dlat);
    double dVt = 1.5 * std::sqrt(3.) * sqr(sech(dRho)) * std::tanh(dRho);
    double dOmega;
    if (dRho == 0.) {
        dOmega = 0.;
    }
    else {
        dOmega = dVt / dRho;
    }

    return 2. * (1. + std::tanh(dRho / dD * std::sin(dlon - dOmega * dT)));
}

double MDPI_gulfstream(double lon, double lat) {
    constexpr double d2r = Constants::degreesToRadians();

    auto sqr                  = [](const double x) { return x * x; };
    auto sech                 = [](const double x) { return 1. / std::cosh(x); };

    constexpr double length     =   1.2 * M_PI;
    constexpr double gf_coef    =   1.0;       // Coefficient for a Gult Stream term (0.0 = no Gulf Stream)
    constexpr double gf_ori_lon = -80.0 * d2r; // Origin of the Gulf Stream (longitude in deg)
    constexpr double gf_ori_lat =  25.0 * d2r; // Origin of the Gulf Stream (latitude in deg)
    constexpr double gf_end_lon =  -1.8 * d2r; // End of the Gulf Stream (longitude in deg)
    constexpr double gf_end_lat =  50.0 * d2r; // End of the Gulf Stream (latitude in deg)
    constexpr double gf_dmp_lon = -25.5 * d2r; // Point of the Gulf Stream decrease (longitude in deg)
    constexpr double gf_dmp_lat = -55.5 * d2r; // Point of the Gulf Stream decrease (latitude in deg)

    double dr0 = std::sqrt(sqr(gf_end_lon - gf_ori_lon) + sqr(gf_end_lat - gf_ori_lat));
    double dr1 = std::sqrt(sqr(gf_dmp_lon - gf_ori_lon) + sqr(gf_dmp_lat - gf_ori_lat));

    double gf_per_lon = [lon,d2r]() {
        double gf_per_lon = lon - 180.;
        while (gf_per_lon > 180.) {
            gf_per_lon -= 360.;
        }
        while (gf_per_lon < -180.) {
            gf_per_lon += 360.;
        }
        return gf_per_lon * d2r;
    }();
    double dx = gf_per_lon - gf_ori_lon;
    double dy = lat * d2r - gf_ori_lat;
    double dr = std::sqrt(sqr(dx) + sqr(dy));
    double dth = std::atan2(dy, dx);
    double dc = 1.3 * gf_coef;
    if (dr > dr0) {
        dc = 0.;
    }
    if (dr > dr1) {
        dc *= std::cos(M_PI_2 * (dr - dr1)/(dr0 - dr1));
    }

    double background_func = MDPI_sinusoid(lon, lat);
    return background_func + dc * (std::max(1000. * std::sin(0.4 * (0.5 * dr + dth) + 0.007 * std::cos(50. * dth) + 0.37 * M_PI), 999.) - 999.);
}

extern "C" {
    double atlas__functions__MDPI_sinusoid(double& lon, double& lat) {
        return MDPI_sinusoid(lon, lat);
    }
    double atlas__functions__MDPI_harmonic(double& lon, double& lat) {
        return MDPI_harmonic(lon, lat);
    }
    double atlas__functions__MDPI_vortex(double& lon, double& lat) {
        return MDPI_vortex(lon, lat);
    }
    double atlas__functions__MDPI_gulfstream(double& lon, double& lat) {
        return MDPI_gulfstream(lon, lat);
    }
}

}  // namespace function
}  // namespace util
}  // namespace atlas
