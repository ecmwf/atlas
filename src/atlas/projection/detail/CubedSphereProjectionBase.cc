/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereProjectionBase.h"

#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>

#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------
// Helper functions and variables local to this translation unit
namespace {

static constexpr double deg2rad = util::Constants::degreesToRadians();

static void schmidtTransform( double stretchFac, double targetLon,
                              double targetLat, double lonlat[]) {


    double c2p1 = 1.0 + stretchFac*stretchFac;
    double c2m1 = 1.0 - stretchFac*stretchFac;

    double sin_p = std::sin(targetLat * deg2rad);
    double cos_p = std::cos(targetLat * deg2rad);

    double sin_lat;
    double cos_lat;
    double lat_t;

    if ( std::abs( c2m1 ) > 1.0e-7 ) {
        sin_lat = std::sin( lonlat[LAT] * deg2rad);
        lat_t   = std::asin( ( c2m1 + c2p1 * sin_lat ) / ( c2p1 + c2m1 * sin_lat ) );
    }
    else {  // no stretching
        lat_t = lonlat[LAT];
    }

    sin_lat      = std::sin( lat_t);
    cos_lat      = std::cos( lat_t);
    double sin_o = -( sin_p * sin_lat + cos_p * cos_lat * cos( lonlat[LON] *deg2rad ) );

    if ( ( 1. - std::abs( sin_o ) ) < 1.0e-7 ) {  // poles
        lonlat[LON] = 0.0;
        lonlat[LAT] = std::copysign( 90.0, sin_o );
    }
    else {
        lonlat[LAT] = std::asin( sin_o );
        lonlat[LON] = targetLon + atan2( -cos_lat * std::sin( lonlat[LON] * deg2rad ),
                                         -sin_lat * cos_p + cos_lat * sin_p * std::cos( lonlat[LON] * deg2rad ) );
        if ( lonlat[LON] < 0.0 ) {
            lonlat[LON] += 360.;
        }
        else if ( lonlat[LON] >= 360. ) {
            lonlat[LON] -= 360.;
        }
    }
}

void sphericalToCartesian( const double lonlat[], double xyz[] ) {
    auto crd_sys            = ProjectionUtilities::CoordinateSystem::LEFT_HAND;
    constexpr double radius = 1.;
    ProjectionUtilities::sphericalToCartesian( lonlat, xyz, crd_sys, radius );
}

void cartesianToSpherical( const double xyz[], double lonlat[] ) {
    auto crd_sys            = ProjectionUtilities::CoordinateSystem::LEFT_HAND;
    constexpr double radius = 0.;  // --> equivalent to radius = norm(xyz)
    ProjectionUtilities::cartesianToSpherical( xyz, lonlat, crd_sys, radius );
}

}  // namespace


// -------------------------------------------------------------------------------------------------

CubedSphereProjectionBase::CubedSphereProjectionBase( const eckit::Parametrisation& params ) :
    CubedSphereTiles_(params)
{
  ATLAS_TRACE( "CubedSphereProjectionBase::CubedSphereProjectionBase" );

    // Shift projection by a longitude
    shiftLon_ = 0.0;
    if ( params.has( "ShiftLon" ) ) {
        params.get( "ShiftLon", shiftLon_ );
        ATLAS_ASSERT( shiftLon_ <= 90.0, "ShiftLon should be <= 90.0 degrees" );
        ATLAS_ASSERT( shiftLon_ >= -90.0, "ShiftLon should be >= -90.0 degrees" );
    }

    // Apply a Schmidt transform
    doSchmidt_  = false;
    stretchFac_ = 0.0;
    targetLon_  = 0.0;
    targetLat_  = 0.0;
    if ( params.has( "DoSchmidt" ) ) {
        params.get( "DoSchmidt", doSchmidt_ );
        if ( doSchmidt_ ) {
            params.get( "StretchFac", stretchFac_ );
            params.get( "TargetLon", targetLon_ );
            params.get( "TargetLat", targetLat_ );
        }
    }
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::hash( eckit::Hash& h ) const {
    // Add stretching options to hash
    h.add( shiftLon_ );
    h.add( doSchmidt_ );
    h.add( stretchFac_ );
    h.add( targetLon_ );
    h.add( targetLat_ );
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat_post( double xyz[], const idx_t& t, double crd[] ) const {
    cartesianToSpherical( xyz, crd );

    if ( crd[LON] < 0. ) {
      crd[LON] += 360.;
    }
    crd[LON] -= 180.;

    // Convert to cartesian
    sphericalToCartesian( crd, xyz );

    // Perform tile specific rotation
    tileRotate(t, xyz);

    // Back to lonlat
    cartesianToSpherical( xyz, crd );

    // Shift longitude
    if ( shiftLon_ != 0.0 ) {
        crd[LON] = crd[LON] + shiftLon_;
        if ( crd[LON] < -180.) {
            crd[LON] = 360. + crd[LON];
        }
        if ( crd[LON] > 180. ) {
            crd[LON] = - 360. + crd[LON];
        }
    }

    // To 0, 360
    if ( crd[LON] < 0.0 ) {
        crd[LON] += 360.;
    }

    // Schmidt transform
    if ( doSchmidt_ ) {
        schmidtTransform( stretchFac_, targetLon_, targetLat_, crd );
    }

    // longitude does not make sense at the poles - set to 0.
    if ( std::abs( std::abs( crd[LAT] ) - 90. ) < 1e-15 ) {
        crd[LON] = 0.;
    }

}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::lonlat2xy_pre( double crd[], idx_t& t, double xyz[] ) const {
    if ( std::abs( crd[LON] ) < 1e-15 ) {
        crd[LON] = 0.;
    }
    if ( std::abs( crd[LAT] ) < 1e-15 ) {
        crd[LAT] = 0.;
    }

    // To [-45.0, 315)
    if ( crd[LON] >= 315.0 ) {
        crd[LON] -= 360.;
    }

    // find tile which this lonlat is linked to
    // works [-45, 315.0)
    t = CubedSphereTiles_.tileFromLonLat(crd);

    sphericalToCartesian(crd, xyz);
    tileRotateInverse(t, xyz);

}

// -------------------------------------------------------------------------------------------------


void CubedSphereProjectionBase::xy2alphabetat( const double xy[], idx_t& t, double ab[] ) const {
    // xy is in degrees while ab is in degree
    // ab are the  (alpha, beta) coordinates and t is the tile index.

    t = CubedSphereTiles_.tileFromXY(xy);
    double normalisedX = xy[XX]/90.;
    double normalisedY = (xy[YY] + 135.)/90.;
    ab[LON] = (normalisedX - CubedSphereTiles_.xy2abOffsets()[XX][t])* 90.0 - 45.0;
    ab[LAT] = (normalisedY - CubedSphereTiles_.xy2abOffsets()[YY][t])* 90.0 - 45.0;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::alphabetat2xy( const idx_t& t, const double ab[], double xy[] ) const {
    // xy and ab are in degrees
    // (alpha, beta) and tiles.

    xy[XX] = ab[LON] + 45.0 + CubedSphereTiles_.ab2xyOffsets()[LON][t];
    xy[YY] = ab[LAT] + 45.0 + CubedSphereTiles_.ab2xyOffsets()[LAT][t];

    CubedSphereTiles_.enforceXYdomain(xy);
}

void CubedSphereProjectionBase::tileRotate(const idx_t& t, double xyz[]) const  {

    switch (t) {
    case 0: CubedSphereTiles_.tile0Rotate(xyz); break;
    case 1: CubedSphereTiles_.tile1Rotate(xyz); break;
    case 2: CubedSphereTiles_.tile2Rotate(xyz); break;
    case 3: CubedSphereTiles_.tile3Rotate(xyz); break;
    case 4: CubedSphereTiles_.tile4Rotate(xyz); break;
    case 5: CubedSphereTiles_.tile5Rotate(xyz); break;
    }
}

void CubedSphereProjectionBase::tileRotateInverse(const idx_t& t, double xyz[]) const  {

    switch (t) {
    case 0: CubedSphereTiles_.tile0RotateInverse(xyz); break;
    case 1: CubedSphereTiles_.tile1RotateInverse(xyz); break;
    case 2: CubedSphereTiles_.tile2RotateInverse(xyz); break;
    case 3: CubedSphereTiles_.tile3RotateInverse(xyz); break;
    case 4: CubedSphereTiles_.tile4RotateInverse(xyz); break;
    case 5: CubedSphereTiles_.tile5RotateInverse(xyz); break;
    }
}

// -------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
