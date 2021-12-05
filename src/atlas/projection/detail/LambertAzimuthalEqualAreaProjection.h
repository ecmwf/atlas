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

class LambertAzimuthalEqualAreaProjection final : public ProjectionImpl {
public:
    // constructor
    LambertAzimuthalEqualAreaProjection(const eckit::Parametrisation&);

    // projection name
    static std::string static_type() { return "lambert_azimuthal_equal_area"; }
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
    PointLonLat reference_;  ///< central longitude/standard parallel [degree]/[degree]
    double radius_;          ///< sphere radius

    double lambda0_;  ///< central longitude [rad]
    double phi1_;     ///< standard parallel [rad]
    double sin_phi1_;
    double cos_phi1_;
    double false_northing_{0};
    double false_easting_{0};
};

}  // namespace detail
}  // namespace projection
}  // namespace atlas
