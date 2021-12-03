/*
 * (C) Copyright 2013 ECMWF.
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
#include "atlas/util/NormaliseLongitude.h"

namespace atlas {
namespace projection {
namespace detail {

template <typename Rotation>
class MercatorProjectionT final : public ProjectionImpl {
public:
    // constructor
    MercatorProjectionT(const eckit::Parametrisation&);

    // projection name
    static std::string static_type() { return Rotation::typePrefix() + "mercator"; }
    std::string type() const override { return static_type(); }

    // projection and inverse projection
    void xy2lonlat(double crd[]) const override;
    void lonlat2xy(double crd[]) const override;

    Jacobian jacobian(const PointLonLat&) const override;

    bool strictlyRegional() const override { return true; }  // Mercator projection cannot be used for global grids
    RectangularLonLatDomain lonlatBoundingBox(const Domain& domain) const override {
        return ProjectionImpl::lonlatBoundingBox(domain);
    }

    // specification
    Spec spec() const override;

    std::string units() const override { return "meters"; }

    void hash(eckit::Hash&) const override;

protected:
    Normalise normalise_;
    util::NormaliseLongitude normalise_mercator_;
    double lon0_;          // central longitude (default = 0 )
    double lat1_;          // latitude where cylinder cuts sphere (default = 0 )
    double radius_;        // sphere radius
    double k_radius_;      // sphere radius
    double inv_k_radius_;  // 1/(sphere radius)

    double eccentricity_;
    double semi_major_axis_;
    double semi_minor_axis_;

    double false_easting_;
    double false_northing_;


    void setup(const eckit::Parametrisation& p);

private:
    Rotation rotation_;
};

using MercatorProjection        = MercatorProjectionT<NotRotated>;
using RotatedMercatorProjection = MercatorProjectionT<Rotated>;

}  // namespace detail
}  // namespace projection
}  // namespace atlas
