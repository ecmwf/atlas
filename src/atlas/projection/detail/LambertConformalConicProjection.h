/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/domain.h"
#include "atlas/projection/detail/ProjectionImpl.h"

namespace atlas {
namespace projection {
namespace detail {

class LambertConformalConicProjection final : public ProjectionImpl {
public:
    // constructor
    LambertConformalConicProjection(const eckit::Parametrisation&);

    // projection name
    static std::string static_type() { return "lambert_conformal_conic"; }
    std::string type() const override { return static_type(); }

    // projection and inverse projection
    void xy2lonlat(double crd[]) const override;
    void lonlat2xy(double crd[]) const override;

    Jacobian jacobian(const PointLonLat&) const override;

    bool strictlyRegional() const override { return true; }

    RectangularLonLatDomain lonlatBoundingBox(const Domain& domain) const override {
        return ProjectionImpl::lonlatBoundingBox(domain);
    }

    // specification
    Spec spec() const override;

    std::string units() const override { return "meters"; }

    void hash(eckit::Hash&) const override;

private:
    double radius_;  ///< sphere radius
    double lat1_;    ///< first latitude from the pole at which the secant cone cuts the sphere
    double lat2_;    ///< second latitude from the pole at which the secant cone cuts the sphere
    double lat0_;    ///< latitude of origin, where Dx and Dy are specified
    double lon0_;    ///< longitude of origin, meridian parallel to y-axis along which latitude increases
                     ///  as the y-coordinate increases

    double F_;
    double n_;
    double inv_n_;
    double rho0_;
    double sign_;
};

}  // namespace detail
}  // namespace projection
}  // namespace atlas
