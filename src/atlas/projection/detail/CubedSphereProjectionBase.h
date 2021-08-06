/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/grid/Tiles.h"
#include "atlas/library/config.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/utils/Hash.h"

namespace atlas {
class CubedSphereTiles;
}

namespace atlas {
namespace projection {
namespace detail {

class CubedSphereProjectionBase : public ProjectionImpl {
public:
    // constructor
    CubedSphereProjectionBase( const eckit::Parametrisation& );

    void hash( eckit::Hash& ) const;

    atlas::grid::CubedSphereTiles getCubedSphereTiles() const { return tiles_; };

protected:
    // projection and inverse projection
    void xy2lonlat_post( double xyz[], const idx_t t, double crd[] ) const;
    void lonlat2xy_pre( double crd[], idx_t& t, double xyz[] ) const;

    void xy2alphabetat( const double xy[], idx_t& t, double ab[] ) const;
    void alphabetat2xy( const idx_t t, const double ab[], double xy[] ) const;

private:
    atlas::grid::CubedSphereTiles tiles_;
    // Shift entire grid
    double shiftLon_;
    // Schmidt transform
    bool doSchmidt_;
    double stretchFac_;
    double targetLon_;
    double targetLat_;

    std::array<std::array<double, 6>, 2> tiles_offsets_ab2xy_;
    std::array<std::array<double, 6>, 2> tiles_offsets_xy2ab_;
};

}  // namespace detail
}  // namespace projection
}  // namespace atlas
