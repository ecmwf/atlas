/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/array.h"
#include "atlas/domain.h"
#include "atlas/projection/detail/CubedSphereProjectionBase.h"
#include "atlas/projection/detail/ProjectionImpl.h"

namespace atlas {
namespace projection {
namespace detail {

class CubedSphereEquiAnglProjection final : public ProjectionImpl, public CubedSphereProjectionBase {
  public:
    // constructor
    CubedSphereEquiAnglProjection( const eckit::Parametrisation& );

    // projection name
    static std::string static_type() { return "cubedsphere_equiangular"; }
    std::string type() const override { return static_type(); }

    // projection and inverse projection
    void xy2lonlat( double crd[] ) const override;
    void lonlat2xy( double crd[] ) const override;

    Jacobian jacobian( const PointLonLat& ) const override;

    bool strictlyRegional() const override { return false; }
    RectangularLonLatDomain lonlatBoundingBox( const Domain& domain ) const override {
        return ProjectionImpl::lonlatBoundingBox( domain );
    }

    // specification
    Spec spec() const override;

    std::string units() const override { return "degrees"; }

    void hash( eckit::Hash& ) const override;
};

}  // namespace detail
}  // namespace projection
}  // namespace atlas
