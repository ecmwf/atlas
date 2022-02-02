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

namespace atlas {
namespace projection {
namespace detail {

template <typename Rotation>
class SchmidtProjectionT final : public ProjectionImpl {
public:
    // constructor
    SchmidtProjectionT(const eckit::Parametrisation& p);
    SchmidtProjectionT();

    // projection name
    static std::string static_type() { return Rotation::typePrefix() + "schmidt"; }
    std::string type() const override { return static_type(); }

    // projection and inverse projection
    void xy2lonlat(double crd[]) const override;
    void lonlat2xy(double crd[]) const override;

    Jacobian jacobian(const PointLonLat&) const override;

    bool strictlyRegional() const override { return false; }  // schmidt is global grid
    RectangularLonLatDomain lonlatBoundingBox(const Domain& domain) const override {
        return ProjectionImpl::lonlatBoundingBox(domain);
    }

    // specification
    Spec spec() const override;

    std::string units() const override { return "degrees"; }

    void hash(eckit::Hash&) const override;

private:
    double c_;  // stretching factor
    Rotation rotation_;
    PointXYZ north0_;
    PointXYZ north1_;
};

using SchmidtProjection        = SchmidtProjectionT<NotRotated>;
using RotatedSchmidtProjection = SchmidtProjectionT<Rotated>;

}  // namespace detail
}  // namespace projection
}  // namespace atlas
