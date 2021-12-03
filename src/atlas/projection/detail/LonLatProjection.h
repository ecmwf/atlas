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

#include "atlas/projection/detail/ProjectionImpl.h"

namespace atlas {
class Projection;
}

namespace atlas {
namespace projection {
namespace detail {

template <typename Rotation>
class LonLatProjectionT final : public ProjectionImpl {
private:
    friend class atlas::Projection;
    LonLatProjectionT() = default;

public:
    // constructor
    LonLatProjectionT(const eckit::Parametrisation&);

    // destructor
    ~LonLatProjectionT() = default;

    // projection name
    static std::string static_type() { return Rotation::typePrefix() + "lonlat"; }
    std::string type() const override { return static_type(); }

    // projection and inverse projection
    void xy2lonlat(double crd[]) const override { rotation_.rotate(crd); }
    void lonlat2xy(double crd[]) const override { rotation_.unrotate(crd); }

    Jacobian jacobian(const PointLonLat&) const override;

    bool strictlyRegional() const override { return false; }
    RectangularLonLatDomain lonlatBoundingBox(const Domain&) const override;

    // specification
    Spec spec() const override;

    std::string units() const override { return "degrees"; }

    operator bool() const override { return rotation_.rotated(); }

    void hash(eckit::Hash&) const override;

private:
    Rotation rotation_;
};

using LonLatProjection        = LonLatProjectionT<NotRotated>;
using RotatedLonLatProjection = LonLatProjectionT<Rotated>;

// --------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
