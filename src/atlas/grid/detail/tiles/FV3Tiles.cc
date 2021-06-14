/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <array>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <iomanip>

#include "eckit/utils/Hash.h"

#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/grid/detail/tiles/FV3Tiles.h"
#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"


namespace {

static constexpr bool debug = false; // constexpr so compiler can optimize `if ( debug ) { ... }` out

static constexpr double deg2rad = ::atlas::util::Constants::degreesToRadians();
static constexpr double rad2deg = ::atlas::util::Constants::radiansToDegrees();

using atlas::projection::detail::ProjectionUtilities;

static bool is_tiny( const double& x ) {
    constexpr double epsilon = 1.e-15;
    return (std::abs(x) < epsilon );
}

static bool is_same( const double& x, const double& y, const double& tol = 1.0 ) {
    constexpr double epsilon = 1.e-15;
    return (std::abs(x-y) < epsilon * tol);
}

void sphericalToCartesian(const double lonlat[], double xyz[] ) {
    auto crd_sys = ProjectionUtilities::CoordinateSystem::LEFT_HAND;
    constexpr double radius = 1.;
    ProjectionUtilities::sphericalToCartesian(lonlat, xyz, crd_sys, radius);
}



}

namespace atlas {
namespace cubedspheretiles {

// constructor
FV3CubedSphereTiles::FV3CubedSphereTiles( const eckit::Parametrisation& ) {

}

idx_t FV3CubedSphereTiles::tileFromXY( const double xy[] ) const  {

    // Assume one face-edge is of length 90 degrees.
    //
    //   y ^
    //     |
    //    135              ----------
    //     |              |     ^    |
    //     |              |          |
    //     |              |=<   2   <|
    //     |              |     v    |
    //     |              |     =    |
    //     45  0----------2----------3----------4----------
    //     |   |    ^     |     ^    |    =     |     =    |
    //     |   |          |          |    ^     |     ^    |
    //     |   |=<  0    <|=<   1   <|=<  3    <|=<   4   <|
    //     |   |    v     |     v    |          |          |
    //     |   |    =     |     =    |    v     |     v    |
    //    -45  0 ---------1----------1----------5----------
    //     |                                    |     =    |
    //     |                                    |     ^    |
    //     |                                    |=<   5   <|
    //     |                                    |          |
    //     |                                    |     v    |
    //   -135                                    ----------(5 for end iterator)
    //     ----0---------90--------180--------270--------360--->  x



    idx_t t{-1};

    if ((xy[XX] >= 0.) && ( xy[YY] >= -45.) && (xy[XX] < 90.) && (xy[YY] < 45.)) {
       t = 0;
    } else if ((xy[XX] >= 90.) && ( xy[YY] >= -45.) && (xy[XX] < 180.) && (xy[YY] < 45.)) {
       t = 1;
    } else if ((xy[XX] >= 90.) && ( xy[YY] >= 45.) && (xy[XX] < 180.) && (xy[YY] < 135.)) {
       t = 2;
    } else if ((xy[XX] >= 180.) && ( xy[YY] > -45.) && (xy[XX] < 270.) && (xy[YY] <= 45.)) {
       t = 3;
    } else if ((xy[XX] >= 270.) && ( xy[YY] > -45.) && (xy[XX] < 360.) && (xy[YY] <= 45.)) {
       t = 4;
    } else if ((xy[XX] >= 270.) && ( xy[YY] > -135.) && (xy[XX] < 360.) && (xy[YY] <= -45.)) {
       t = 5;
    }

    // extra points
    if ( is_same( xy[XX], 0.   ) && is_same( xy[YY],  45. ) ) { t = 0; }
    if ( is_same( xy[XX], 180. ) && is_same( xy[YY], -45. ) ) { t = 1; }

    // for end iterator !!!!
    if ( is_same( xy[XX], 360. ) && is_same( xy[YY], -135.) ) { t = 5; }

    ATLAS_ASSERT( t >= 0 );

    return t;
}

idx_t FV3CubedSphereTiles::tileFromLonLat( const double crd[] ) const {

    idx_t t(-1); // tile index
    double xyz[3];

    sphericalToCartesian(crd, xyz);

    static const double cornerLat= std::asin(1./std::sqrt(3.0));
      // the magnitude of the latitude at the corners of the cube (not including the sign)
      // in radians.

    const double & lon = crd[LON];
    const double & lat = crd[LAT];

    double zPlusAbsX = xyz[ZZ] + std::abs(xyz[XX]);
    double zPlusAbsY = xyz[ZZ] + std::abs(xyz[YY]);
    double zMinusAbsX = xyz[ZZ] - std::abs(xyz[XX]);
    double zMinusAbsY = xyz[ZZ] - std::abs(xyz[YY]);

    // Note that this method can lead to roundoff errors that can
    // cause the tile selection to fail.
    // To this end we enforce that tiny values close (in roundoff terms)
    // to a boundary should end up exactly on the boundary.
    if ( is_tiny(zPlusAbsX) ) { zPlusAbsX = 0.; }
    if ( is_tiny(zPlusAbsY) ) { zPlusAbsY = 0.; }
    if ( is_tiny(zMinusAbsX) ) { zMinusAbsX = 0.; }
    if ( is_tiny(zMinusAbsY) ) { zMinusAbsY = 0.; }

    if (lon >= 1.75 * M_PI  || lon < 0.25 * M_PI) {
        if  ( (zPlusAbsX <= 0.) && (zPlusAbsY <= 0.) ) {
           t = 2;
        } else if ( (zMinusAbsX > 0.) && (zMinusAbsY > 0.) ) {
           t = 5;
        } else {
           t = 0;
        }
        // extra point corner point
        if ( is_same( lon, -0.25 * M_PI ) && is_same( lat, cornerLat ) ) { t = 0; }
    }

    if (lon >= 0.25 * M_PI  && lon < 0.75 * M_PI) {
        // interior
        if  ( (zPlusAbsX <= 0.) && (zPlusAbsY <= 0.) ) {
            t = 2;
        } else if  ( (zMinusAbsX > 0.) && (zMinusAbsY > 0.) ) {
            t = 5;
        } else {
            t = 1;
        }
    }

    if (lon >= 0.75 * M_PI  && lon < 1.25 * M_PI) {
        // interior
        if  ( (zPlusAbsX < 0.) && (zPlusAbsY < 0.) ) {
            t = 2;
        } else if  ( (zMinusAbsX >= 0.) && (zMinusAbsY >= 0.) ) {
            t = 5;
        } else {
            t = 3;
        }
        // extra point corner point
        if ( is_same(lon, 0.75 * M_PI) && is_same( lat, -cornerLat ) ) { t = 1; }
    }

    if (lon >= 1.25 * M_PI  && lon < 1.75 * M_PI) {
        // interior
        if  ( (zPlusAbsX < 0.) && (zPlusAbsY < 0.) ) {
            t = 2;
        } else if  ( (zMinusAbsX >= 0.) && (zMinusAbsY >= 0.) ) {
            t = 5;
        } else {
            t = 4;
        }
    }

    if( debug ) {
        Log::info() << "tileFromLonLat:: lonlat abs xyz t = "
                     << std::setprecision(std::numeric_limits<double>::digits10 + 1)
                     << zPlusAbsX << " " << zPlusAbsY << " " << zMinusAbsX << " " << zMinusAbsY << " "
                     << crd[LON] << " " << crd[LAT] << " "
                     << xyz[XX] << " " << xyz[YY] << " " << xyz[ZZ] << " " << t
                     << std::endl;
    }

    return t;
}

void FV3CubedSphereTiles::enforceXYdomain( double xy[] ) const {
    // the conversion from lonlat to xy can involve some small errors and small errors
    // can affect whether xy that is created within a valid space
    // This has been tested with N = 512 with equiangular and equidistant projections.
    const double tol{70.0};
    constexpr double epsilon = std::numeric_limits<double>::epsilon();


    if ( debug ) {
        Log::info() << "enforcXYDomain before " << xy[XX] << " " << xy[YY] << std::endl;
    }

    xy[XX] = std::max(xy[XX], 0.0);
    xy[XX] = std::min(xy[XX], 360.0 - epsilon);
    xy[YY] = std::max(xy[YY], -135.0 + epsilon);
    xy[YY] = std::min(xy[YY], 135.0 - epsilon);
    if (is_same(xy[XX], 90.0, tol)) { xy[XX] = 90.0; }
    if (is_same(xy[XX], 180.0, tol)) { xy[XX] = 180.0; }
    if (is_same(xy[XX], 270.0, tol)) { xy[XX] = 270.0; }
    if (is_same(xy[YY], -45.0, tol) && (xy[XX] <= 180.0)) { xy[YY] = -45.0; }
    if (is_same(xy[YY], 45.0, tol) && (xy[XX] >= 180.0)) { xy[YY] = 45.0; }
    if (is_same(xy[YY], 45.0, tol) && is_same(xy[XX], 0.0, tol)) {
        xy[XX] = 0.0;
        xy[YY] = 45.0;
    }

    if ( debug ) {
        Log::info() << "enforcXYDomain after " << xy[XX] << " " << xy[YY] << std::endl;
    }
}

FV3CubedSphereTiles::Spec FV3CubedSphereTiles::spec() const {
    Spec tile_spec;
    return tile_spec;
}

void FV3CubedSphereTiles::print( std::ostream& os) const {
    os << "CubedSphereTiles["
       << "]";
}

void FV3CubedSphereTiles::hash( eckit::Hash& ) const {

}

}  // namespace cubedspheretiles
}  // namespace atlas
