/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/domain.h"
#include "atlas/projection/detail/CubedSphereProjectionBase.h"
#include "atlas/projection/detail/ProjectionImpl.h"

namespace atlas {
namespace projection {
namespace detail {

class CubedSphereEquiDistProjection final : public CubedSphereProjectionBase {
public:
    // constructor
    CubedSphereEquiDistProjection(const eckit::Parametrisation&);

    virtual ~CubedSphereEquiDistProjection() {}

    // projection name
    static std::string static_type() { return "cubedsphere_equidistant"; }
    std::string type() const override { return static_type(); }

    /// @brief Convert (x, y) coordinate to (alpha, beta) on tile t.
    void xy2alphabeta(double crd[], idx_t t) const override;

    /// @brief Convert (alpha, beta) coordinate to (x, y) on tile t.
    void alphabeta2xy(double crd[], idx_t t) const override;

    /// @brief Convert (lon, lat) coordinate to (alpha, beta) on tile t.
    void lonlat2alphabeta(double crd[], idx_t t) const override;

    /// @brief Convert (alpha, beta) coordinate to (lon, lat) on tile t.
    void alphabeta2lonlat(double crd[], idx_t t) const override;

    /// @brief Jacobian of (alpha, beta) with respect to (lon, lat) on tile t
    Jacobian alphabetaJacobian(const PointLonLat& lonlat, idx_t t) const override;

    /// @brief Jacobian of (x, y) with respect to (lon, lat) on tile t
    Jacobian jacobian(const PointLonLat& lonlat, idx_t t) const override;

    // projection and inverse projection
    void xy2lonlat(double crd[]) const override;
    void lonlat2xy(double crd[]) const override;

    Jacobian jacobian(const PointLonLat&) const override;

    bool strictlyRegional() const override { return false; }
    RectangularLonLatDomain lonlatBoundingBox(const Domain& domain) const override {
        return ProjectionImpl::lonlatBoundingBox(domain);
    }

    // specification
    Spec spec() const override;

    std::string units() const override { return "degrees"; }

    void hash(eckit::Hash&) const override;
};

}  // namespace detail
}  // namespace projection
}  // namespace atlas
