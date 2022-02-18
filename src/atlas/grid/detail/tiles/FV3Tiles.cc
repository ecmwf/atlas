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
#include <iomanip>
#include <limits>
#include <ostream>
#include <string>
#include <utility>

#include "eckit/utils/Hash.h"

#include "atlas/grid/detail/tiles/TilesFactory.h"

#include "atlas/grid/detail/tiles/FV3Tiles.h"
#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"


namespace {

static constexpr bool debug     = false;  // constexpr so compiler can optimize `if ( debug ) { ... }` out
static constexpr double rad2deg = atlas::util::Constants::radiansToDegrees();

using atlas::projection::detail::ProjectionUtilities;

static bool is_tiny(const double& x) {
    constexpr double epsilon = 1.e-12;
    return (std::abs(x) < epsilon);
}

static bool is_same(const double& x, const double& y, const double& tol = 1.0) {
    constexpr double epsilon = 1.e-12;
    return (std::abs(x - y) < epsilon * tol);
}

static bool is_less(const double& lhs, const double& rhs) {
    constexpr double epsilon = 1.e-12;
    return lhs < rhs - epsilon;
}

static bool is_geq(const double& lhs, const double& rhs) {
    constexpr double epsilon = 1.e-12;
    return lhs >= rhs - epsilon;
}

void sphericalToCartesian(const double lonlat[], double xyz[]) {
    auto crd_sys            = ProjectionUtilities::CoordinateSystem::LEFT_HAND;
    constexpr double radius = 1.;
    ProjectionUtilities::sphericalToCartesian(lonlat, xyz, crd_sys, radius);
}

}  // namespace

namespace atlas {
namespace grid {
namespace detail {

// constructor
FV3CubedSphereTiles::FV3CubedSphereTiles(const eckit::Parametrisation&) {}

std::array<std::array<double, 6>, 2> FV3CubedSphereTiles::xy2abOffsets() const {
    return {{{0., 1., 1., 2., 3., 3.}, {1., 1., 2., 1., 1., 0.}}};
}

std::array<std::array<double, 6>, 2> FV3CubedSphereTiles::ab2xyOffsets() const {
    return {{{0., 90., 90., 180., 270., 270.}, {-45., -45., 45., -45., -45., -135.}}};
}

static void tile0Rotate(double xyz[]) {
    //  Face 0, no rotation.
}

static void tile1Rotate(double xyz[]) {
    //  Face 1: rotate -90.0 degrees about z axis
    constexpr double angle = -M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}

static void tile2Rotate(double xyz[]) {
    //  Face 2: rotate -90.0 degrees about z axis
    //          rotate  90.0 degrees about x axis
    constexpr double angle_z = -M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle_z, xyz);
    constexpr double angle_x = M_PI / 2.0;
    ProjectionUtilities::rotate3dX(angle_x, xyz);
}

static void tile3Rotate(double xyz[]) {
    //  Face 3: rotate -180.0 degrees about z axis
    constexpr double angle = -M_PI;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}

static void tile4Rotate(double xyz[]) {
    //  Face 4: rotate 90.0 degrees about z axis
    constexpr double angle = M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}

static void tile5Rotate(double xyz[]) {
    // Face 5: rotate 90.0 degrees about y axis
    //         rotate 90.0 degrees about z axis
    constexpr double angle_y = M_PI / 2.0;
    ProjectionUtilities::rotate3dY(angle_y, xyz);
    constexpr double angle_z = M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle_z, xyz);
}

void FV3CubedSphereTiles::rotate(idx_t t, double xyz[]) const {
    switch (t) {
        case 0: {
            tile0Rotate(xyz);
            break;
        }
        case 1: {
            tile1Rotate(xyz);
            break;
        }
        case 2: {
            tile2Rotate(xyz);
            break;
        }
        case 3: {
            tile3Rotate(xyz);
            break;
        }
        case 4: {
            tile4Rotate(xyz);
            break;
        }
        case 5: {
            tile5Rotate(xyz);
            break;
        }
        default: {
            throw_OutOfRange("t", t, 6);
        }
    }
}

static void tile0RotateInverse(double xyz[]) {
    //  Face 0, no rotation.
}

static void tile1RotateInverse(double xyz[]) {
    //  Face 1: rotate 90.0 degrees about z axis
    constexpr double angle = M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}

static void tile2RotateInverse(double xyz[]) {
    //  Face 2: rotate  90.0 degrees about x axis
    //          rotate -90.0 degrees about z axis
    constexpr double angle_x = -M_PI / 2.0;
    ProjectionUtilities::rotate3dX(angle_x, xyz);
    constexpr double angle_z = M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle_z, xyz);
}

static void tile3RotateInverse(double xyz[]) {
    //  Face 3: rotate  180.0 degrees about z axis
    constexpr double angle = M_PI;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}

static void tile4RotateInverse(double xyz[]) {
    //  Face 4: rotate -90.0 degrees about y axis
    constexpr double angle = -M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}

static void tile5RotateInverse(double xyz[]) {
    //  Face 5: rotate -90.0 degrees about y axis
    //          rotate -90.0 degrees about z axis
    constexpr double angle_z = -M_PI / 2.;
    ProjectionUtilities::rotate3dZ(angle_z, xyz);
    constexpr double angle_y = -M_PI / 2.;
    ProjectionUtilities::rotate3dY(angle_y, xyz);
}

void FV3CubedSphereTiles::unrotate(idx_t t, double xyz[]) const {
    switch (t) {
        case 0: {
            tile0RotateInverse(xyz);
            break;
        }
        case 1: {
            tile1RotateInverse(xyz);
            break;
        }
        case 2: {
            tile2RotateInverse(xyz);
            break;
        }
        case 3: {
            tile3RotateInverse(xyz);
            break;
        }
        case 4: {
            tile4RotateInverse(xyz);
            break;
        }
        case 5: {
            tile5RotateInverse(xyz);
            break;
        }
        default: {
            throw_OutOfRange("t", t, 6);
        }
    }
}

idx_t FV3CubedSphereTiles::indexFromXY(const double xy[]) const {
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

    if ((xy[XX] >= 0.) && (xy[YY] >= -45.) && (xy[XX] < 90.) && (xy[YY] < 45.)) {
        t = 0;
    }
    else if ((xy[XX] >= 90.) && (xy[YY] >= -45.) && (xy[XX] < 180.) && (xy[YY] < 45.)) {
        t = 1;
    }
    else if ((xy[XX] >= 90.) && (xy[YY] >= 45.) && (xy[XX] < 180.) && (xy[YY] < 135.)) {
        t = 2;
    }
    else if ((xy[XX] >= 180.) && (xy[YY] > -45.) && (xy[XX] < 270.) && (xy[YY] <= 45.)) {
        t = 3;
    }
    else if ((xy[XX] >= 270.) && (xy[YY] > -45.) && (xy[XX] < 360.) && (xy[YY] <= 45.)) {
        t = 4;
    }
    else if ((xy[XX] >= 270.) && (xy[YY] > -135.) && (xy[XX] < 360.) && (xy[YY] <= -45.)) {
        t = 5;
    }

    // extra points
    if (is_same(xy[XX], 0.) && is_same(xy[YY], 45.)) {
        t = 0;
    }
    if (is_same(xy[XX], 180.) && is_same(xy[YY], -45.)) {
        t = 1;
    }

    // for end iterator !!!!
    if (is_same(xy[XX], 360.) && is_same(xy[YY], -135.)) {
        t = 5;
    }

    ATLAS_ASSERT(t >= 0);

    return t;
}

// Calculates the FV3 panel
// Input (crd) is Longitude and Latitude in Radians
// Output is tile number
idx_t FV3CubedSphereTiles::indexFromLonLat(const double crd[]) const {
    idx_t t(-1);  // tile index
    double xyz[3];

    sphericalToCartesian(crd, xyz);

    static const double cornerLat = std::asin(1. / std::sqrt(3.0)) * rad2deg;
    // the magnitude of the latitude at the corners of the cube (not including the sign)
    // in radians.

    const double& lon = crd[LON];
    const double& lat = crd[LAT];

    double zPlusAbsX  = xyz[ZZ] + std::abs(xyz[XX]);
    double zPlusAbsY  = xyz[ZZ] + std::abs(xyz[YY]);
    double zMinusAbsX = xyz[ZZ] - std::abs(xyz[XX]);
    double zMinusAbsY = xyz[ZZ] - std::abs(xyz[YY]);

    // Note that this method can lead to roundoff errors that can
    // cause the tile selection to fail.
    // To this end we enforce that tiny values close (in roundoff terms)
    // to a boundary should end up exactly on the boundary.
    if (is_tiny(zPlusAbsX)) {
        zPlusAbsX = 0.;
    }
    if (is_tiny(zPlusAbsY)) {
        zPlusAbsY = 0.;
    }
    if (is_tiny(zMinusAbsX)) {
        zMinusAbsX = 0.;
    }
    if (is_tiny(zMinusAbsY)) {
        zMinusAbsY = 0.;
    }

    if (is_geq(lon, 315.) || is_less(lon, 45.)) {
        if ((zPlusAbsX <= 0.) && (zPlusAbsY <= 0.)) {
            t = 2;
        }
        else if ((zMinusAbsX > 0.) && (zMinusAbsY > 0.)) {
            t = 5;
        }
        else {
            t = 0;
        }
        // extra point corner point (considering two different longitudes depending on
        // whether we are on longitude points [0,360] or [-45.,45.]
        if (is_same(lon, -45.) && is_same(lat, cornerLat)) {
            t = 0;
        }
        if (is_same(lon, 315.) && is_same(lat, cornerLat)) {
            t = 0;
        }
    }

    if (is_geq(lon, 45.) && is_less(lon, 135.)) {
        // interior
        if ((zPlusAbsX <= 0.) && (zPlusAbsY <= 0.)) {
            t = 2;
        }
        else if ((zMinusAbsX > 0.) && (zMinusAbsY > 0.)) {
            t = 5;
        }
        else {
            t = 1;
        }
    }

    if (is_geq(lon, 135.) && is_less(lon, 225.)) {
        // interior
        if ((zPlusAbsX < 0.) && (zPlusAbsY < 0.)) {
            t = 2;
        }
        else if ((zMinusAbsX >= 0.) && (zMinusAbsY >= 0.)) {
            t = 5;
        }
        else {
            t = 3;
        }
        // extra point corner point
        if (is_same(lon, 135.) && is_same(lat, -cornerLat)) {
            t = 1;
        }
    }

    if (is_geq(lon, 225.) && is_less(lon, 315.)) {
        // interior
        if ((zPlusAbsX < 0.) && (zPlusAbsY < 0.)) {
            t = 2;
        }
        else if ((zMinusAbsX >= 0.) && (zMinusAbsY >= 0.)) {
            t = 5;
        }
        else {
            t = 4;
        }
    }

    if (debug) {
        Log::info() << "tileFromLonLat:: lonlat abs xyz t = "
                    << std::setprecision(std::numeric_limits<double>::digits10 + 1) << zPlusAbsX << " " << zPlusAbsY
                    << " " << zMinusAbsX << " " << zMinusAbsY << " " << crd[LON] << " " << crd[LAT] << " " << xyz[XX]
                    << " " << xyz[YY] << " " << xyz[ZZ] << " " << t << std::endl;
    }

    return t;
}

void FV3CubedSphereTiles::enforceXYdomain(double xy[]) const {
    // the conversion from lonlat to xy can involve some small errors and small errors
    // can affect whether xy that is created within a valid space
    // This has been tested with N = 512 with equiangular and equidistant projections.
    const double tol{70.0};
    constexpr double epsilon = std::numeric_limits<double>::epsilon();


    if (debug) {
        Log::info() << "enforceXYDomain before " << xy[XX] << " " << xy[YY] << std::endl;
    }

    xy[XX] = std::max(xy[XX], 0.0);
    xy[XX] = std::min(xy[XX], 360.0 - epsilon);
    xy[YY] = std::max(xy[YY], -135.0 + epsilon);
    xy[YY] = std::min(xy[YY], 135.0 - epsilon);
    if (is_same(xy[XX], 90.0, tol)) {
        xy[XX] = 90.0;
    }
    if (is_same(xy[XX], 180.0, tol)) {
        xy[XX] = 180.0;
    }
    if (is_same(xy[XX], 270.0, tol)) {
        xy[XX] = 270.0;
    }
    if (is_same(xy[YY], -45.0, tol) && (xy[XX] <= 180.0)) {
        xy[YY] = -45.0;
    }
    if (is_same(xy[YY], 45.0, tol) && (xy[XX] >= 180.0)) {
        xy[YY] = 45.0;
    }
    if (is_same(xy[YY], 45.0, tol) && is_same(xy[XX], 0.0, tol)) {
        xy[XX] = 0.0;
        xy[YY] = 45.0;
    }

    if (debug) {
        Log::info() << "enforceXYDomain after " << xy[XX] << " " << xy[YY] << std::endl;
    }
}

atlas::PointXY FV3CubedSphereTiles::tileCubePeriodicity(const atlas::PointXY& xyExtended,
                                                        const atlas::idx_t tile) const {
    throw_NotImplemented("tileCubePeriodicity not implemented for FV3Tiles", Here());
}

void FV3CubedSphereTiles::print(std::ostream& os) const {
    os << "FV3CubedSphereTiles"
       << "]";
}


const PointXY& FV3CubedSphereTiles::tileCentre(size_t t) const {
    throw_NotImplemented("tileCentre not implemented for FV3Tiles", Here());
}

const Jacobian& FV3CubedSphereTiles::tileJacobian(size_t t) const {
    throw_NotImplemented("tileJacobian not implemented for FV3Tiles", Here());;
}


namespace {
static CubedSphereTilesBuilder<FV3CubedSphereTiles> register_builder(FV3CubedSphereTiles::static_type());
}


}  // namespace detail
}  // namespace grid
}  // namespace atlas
