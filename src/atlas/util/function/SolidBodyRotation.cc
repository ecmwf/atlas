/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>

#include "atlas/util/Constants.h"

#include "atlas/util/function/SolidBodyRotation.h"

namespace atlas {

namespace util {

namespace function {

SolidBodyRotation::SolidBodyRotation(const double beta, const double radius) {
    sin_beta_ = std::sin(beta * Constants::degreesToRadians());
    cos_beta_ = std::cos(beta * Constants::degreesToRadians());
    radius_ = radius;
}

void SolidBodyRotation::wind(const double lon, const double lat, double& u, double& v) const {
    double x = lon * Constants::degreesToRadians();
    double y = lat * Constants::degreesToRadians();
    double cos_x = std::cos(x);
    double cos_y = std::cos(y);
    double sin_x = std::sin(x);
    double sin_y = std::sin(y);
    u            = cos_y * cos_beta_ + cos_x * sin_y * sin_beta_;
    v            = -sin_x * sin_beta_;
}

double SolidBodyRotation::u(const double lon, const double lat) const {
    double x = lon * Constants::degreesToRadians();
    double y = lat * Constants::degreesToRadians();
    double cos_x = std::cos(x);
    double cos_y = std::cos(y);
    double sin_y = std::sin(y);
    return cos_y * cos_beta_ + cos_x * sin_y * sin_beta_;
}

double SolidBodyRotation::v(const double lon, const double lat) const {
    double x = lon * Constants::degreesToRadians();
    double sin_x = std::sin(x);
    return -sin_x * sin_beta_;
}

void SolidBodyRotation::vordiv(const double lon, const double lat, double& vor, double& div) const {
    double x = lon * Constants::degreesToRadians();
    double y = lat * Constants::degreesToRadians();

    double cos_x = std::cos(x);
    double cos_y = std::cos(y);
    double sin_x = std::sin(x);
    double sin_y = std::sin(y);

    // Divergence = 1./(R*cos(y)) * ( d/dx( u ) + d/dy( v * cos(y) ) )
    // Vorticity  = 1./(R*cos(y)) * ( d/dx( v ) - d/dy( u * cos(y) ) )
    double ddx_u      = -sin_x * sin_y * sin_beta_;
    double ddy_cosy_v = (-sin_x * sin_beta_) * (-sin_y);
    double ddx_v      = -cos_x * sin_beta_;
    double ddy_cosy_u =
        2 * cos_y * (-sin_y) * cos_beta_ + (-sin_y) * cos_x * sin_y * sin_beta_ + cos_y * cos_x * cos_y * sin_beta_;

    double metric = 1. / (radius_ * cos_y);

    div = metric * (ddx_u + ddy_cosy_v);
    vor = metric * (ddx_v - ddy_cosy_u);
}

double SolidBodyRotation::vorticity(const double lon, const double lat) const {
    double x = lon * Constants::degreesToRadians();
    double y = lat * Constants::degreesToRadians();

    double cos_x = std::cos(x);
    double cos_y = std::cos(y);
    double sin_y = std::sin(y);

    // Vorticity  = 1./(R*cos(y)) * ( d/dx( v ) - d/dy( u * cos(y) ) )
    double ddx_v      = -cos_x * sin_beta_;
    double ddy_cosy_u =
        2 * cos_y * (-sin_y) * cos_beta_ + (-sin_y) * cos_x * sin_y * sin_beta_ + cos_y * cos_x * cos_y * sin_beta_;

    double metric = 1. / (radius_ * cos_y);

    return metric * (ddx_v - ddy_cosy_u);
}

double SolidBodyRotation::divergence(const double lon, const double lat) const {
    double x = lon * Constants::degreesToRadians();
    double y = lat * Constants::degreesToRadians();

    double cos_y = std::cos(y);
    double sin_x = std::sin(x);
    double sin_y = std::sin(y);

    // Divergence = 1./(R*cos(y)) * ( d/dx( u ) + d/dy( v * cos(y) ) )
    double ddx_u      = -sin_x * sin_y * sin_beta_;
    double ddy_cosy_v = (-sin_x * sin_beta_) * (-sin_y);

    double metric = 1. / (radius_ * cos_y);

    return metric * (ddx_u + ddy_cosy_v);
}

double SolidBodyRotation::windMagnitude(const double lon, const double lat) const {
    return std::sqrt(windMagnitudeSquared(lon, lat));
}

double SolidBodyRotation::windMagnitudeSquared(const double lon, const double lat) const {
    double u, v;
    wind(lon, lat, u, v);
    return u * u + v * v;
}

void SolidBodyRotation::windMagnitudeSquaredGradient(const double lon, const double lat, double& dfdx, double& dfdy) const {
    double x = lon * Constants::degreesToRadians();
    double y = lat * Constants::degreesToRadians();

    double cos_x = std::cos(x);
    double cos_y = std::cos(y);
    double sin_x = std::sin(x);
    double sin_y = std::sin(y);

    double metric_y = 1. / radius_;
    double metric_x = metric_y / cos_y;

    double u    = cos_y * cos_beta_ + cos_x * sin_y * sin_beta_;
    double v    = -sin_x * sin_beta_;
    double dudx = metric_x * (-sin_x * sin_y * sin_beta_);
    double dudy = metric_y * (-sin_y * cos_beta_ + cos_x * cos_y * sin_beta_);
    double dvdx = metric_x * (-cos_x * sin_beta_);
    double dvdy = metric_y * (0.);
    dfdx        = 2 * u * dudx + 2 * v * dvdx;
    dfdy        = 2 * u * dudy + 2 * v * dvdy;
}

}  // namespace function

}  // namespace util

}  // namespace atlas
