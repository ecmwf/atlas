/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <ostream>
#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/grid/detail/tiles/LFRicTiles.h"
#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace cubedspheretiles {

namespace {

static constexpr bool debug = false; // constexpr so compiler can optimize `if ( debug ) { ... }` out

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

// constructor
LFRicCubedSphereTiles::LFRicCubedSphereTiles( const eckit::Parametrisation& ) {

}

idx_t LFRicCubedSphereTiles::tileFromXY( const double xy[] ) const  {

    // Assume one face-edge is of length 90 degrees.
    //
    //   y ^
    //     |
    //    135   ----------
    //     |   |    <=    |
    //     |   |          |
    //     |   |=<   4  <=|
    //     |   |     v    |
    //     |   |    <=    |
    //     45  4----------4----------4----------4----------
    //     |   |    ^     |     ^    |    ^     |     ^    |
    //     |   |          |          |          |          |
    //     |   |=<  0    <|=<   1   <|=<  2    <|=<   3   <|
    //     |   |    v     |     v    |    v     |     v    |
    //     |   |    =     |     =    |     =    |     =    |
    //    -45  0 ---------1----------2----------3----------
    //     |   |     ^    |
    //     |   |          |
    //     |   | <   5   <|
    //     |   |          |
    //     |   |     v    |
    //   -135   ----------(5 for end iterator)
    //     ----0---------90--------180--------270--------360--->  x

    idx_t t{-1}; // tile index

    if ((xy[XX] >= 0.) && ( xy[YY] >= -45.) && (xy[XX] < 90.) && (xy[YY] < 45.)) {
       t = 0;
    } else if ((xy[XX] >= 90.) && ( xy[YY] >= -45.) && (xy[XX] < 180.) && (xy[YY] < 45.)) {
       t = 1;
    } else if ((xy[XX] >= 180.) && ( xy[YY] >= -45.) && (xy[XX] < 270.) && (xy[YY] < 45.)) {
       t = 2;
    } else if ((xy[XX] >= 270.) && ( xy[YY] >= -45.) && (xy[XX] < 360.) && (xy[YY] < 45.)) {
       t = 3;
    } else if ((xy[XX] >= 0.) && ( xy[YY] >= 45.) && (xy[XX] <= 90.) && (xy[YY] <= 135.)) {
       t = 4;
    } else if ((xy[XX] > 0.) && ( xy[YY] > -135.) && (xy[XX] < 90.) && (xy[YY] < -45.)) {
       t = 5;
    }

    return t;
}

idx_t LFRicCubedSphereTiles::tileFromLonLat( const double crd[] ) const {
    idx_t t{-1}; // tile index

    double xyz[3];

    sphericalToCartesian(crd, xyz);

    const double & lon = crd[LON];

    double zPlusAbsX = xyz[ZZ] + std::abs(xyz[XX]);
    double zPlusAbsY = xyz[ZZ] + std::abs(xyz[YY]);
    double zMinusAbsX = xyz[ZZ] - std::abs(xyz[XX]);
    double zMinusAbsY = xyz[ZZ] - std::abs(xyz[YY]);

    if ( is_tiny(zPlusAbsX) ) { zPlusAbsX = 0.; }
    if ( is_tiny(zPlusAbsY) ) { zPlusAbsY = 0.; }
    if ( is_tiny(zMinusAbsX) ) { zMinusAbsX = 0.; }
    if ( is_tiny(zMinusAbsY) ) { zMinusAbsY = 0.; }

    if  ( (zPlusAbsX <= 0.) && (zPlusAbsY <= 0.) ) {
         t = 4;
    } else if ( (zMinusAbsX > 0.) && (zMinusAbsY > 0.) ) {
         t = 5;
    } else if (lon >= 1.75 * M_PI  || lon < 0.25 * M_PI) {
         t = 0;
    } else if (lon >= 0.25 * M_PI  && lon < 0.75 * M_PI) {
         t = 1;
    } else if (lon >= 0.75 * M_PI  && lon < 1.25 * M_PI) {
         t = 2;
    } else if (lon >= 1.25 * M_PI  && lon < 1.75 * M_PI) {
         t = 3;
    }

    return t;
}

void LFRicCubedSphereTiles::enforceXYdomain( double xy[] ) const {
    // the conversion from lonlat to xy can involve some small errors and small errors
    // can affect whether xy that is created within a valid space
    // This has been tested with N = 512 with equiangular and equidistant projections.
    const double tol{70.0};
    constexpr double epsilon = std::numeric_limits<double>::epsilon();

    if ( debug ) {
        Log::info() << "enforcXYDomain before " << xy[XX] << " " << xy[YY] << std::endl;
    }

    if (is_same(xy[YY], 45.0, tol)) { xy[YY] = 45.0;}
    if (is_same(xy[YY], -45.0, tol)) { xy[YY] = -45.0;}
    if (xy[YY] < -45.0) {
       xy[XX] = std::min(xy[XX], 90.0 - epsilon);
       xy[XX] = std::max(xy[XX], 0.0 + epsilon);
    } else if (xy[YY] >= -45.0) {
       xy[XX] = std::max(xy[XX], 0.0);
    }
    if (xy[YY] >=  45.0) {
       xy[XX] = std::min(xy[XX], 90.0);
    }

    xy[XX] = std::max(xy[XX], 0.0);
    xy[XX] = std::min(xy[XX], 360.0 - epsilon);
    xy[YY] = std::max(xy[YY], -135.0 + epsilon);
    xy[YY] = std::min(xy[YY], 135.0);


    if ( debug ) {
        Log::info() << "enforcXYDomain after " << xy[XX] << " " << xy[YY] << std::endl;
    }
}

LFRicCubedSphereTiles::Spec LFRicCubedSphereTiles::spec() const {
    Spec tile_spec;
    return tile_spec;
}

void LFRicCubedSphereTiles::print( std::ostream& os) const {
    os << "CubedSphereTiles["
       << "]";
}

void LFRicCubedSphereTiles::hash( eckit::Hash& ) const {

}


}  // namespace cubespheretiles
}  // namespace atlas
