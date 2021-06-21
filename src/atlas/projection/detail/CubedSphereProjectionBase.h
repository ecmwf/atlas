/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once
#include <functional>
#include <memory>

#include "eckit/config/Parametrisation.h"
#include "eckit/utils/Hash.h"
#include "atlas/grid/Tiles.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/library/config.h"

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

    atlas::CubedSphereTiles getCubedSphereTiles() const {return CubedSphereTiles_;};

  protected:

    // projection and inverse projection
    void xy2lonlat_post( double xyz[], const idx_t& t, double crd[] ) const;
    void lonlat2xy_pre( double crd[], idx_t& t, double xyz[] ) const;

    void xy2alphabetat(const double xy[], idx_t & t, double ab[]) const;
    void alphabetat2xy(const idx_t & t, const double ab[],  double xy[]) const;

  private:
    atlas::CubedSphereTiles CubedSphereTiles_;
    // Shift entire grid
    double shiftLon_;
    // Schmidt transform
    bool doSchmidt_;
    double stretchFac_;
    double targetLon_;
    double targetLat_;

    void tileRotate(const idx_t& t, double xyz[]) const;
    void tileRotateInverse(const idx_t& t, double xyz[]) const;

};


}  // namespace detail
}  // namespace projection
}  // namespace atlas
