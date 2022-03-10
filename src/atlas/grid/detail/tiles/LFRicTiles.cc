/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <iostream>
#include <ostream>

#include "atlas/grid/detail/tiles/LFRicTiles.h"
#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/grid/detail/tiles/TilesFactory.h"
#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace grid {
namespace detail {

namespace {

static constexpr bool debug     = false;  // constexpr so compiler can optimize `if ( debug ) { ... }` out
static constexpr double epsilon = 1.e-12;

using projection::detail::ProjectionUtilities;

static bool is_tiny(const double& x) {
    return (std::abs(x) < epsilon);
}

static bool is_same(const double& x, const double& y, const double& tol = 1.0) {
    return (std::abs(x - y) < epsilon * tol);
}

static bool is_less(const double& lhs, const double& rhs) {
    return lhs < rhs - epsilon;
}

static bool is_geq(const double& lhs, const double& rhs) {
    return lhs >= rhs - epsilon;
}

static void sphericalToCartesian(const double lonlat[], double xyz[]) {
    auto crd_sys            = ProjectionUtilities::CoordinateSystem::LEFT_HAND;
    constexpr double radius = 1.;
    ProjectionUtilities::sphericalToCartesian(lonlat, xyz, crd_sys, radius);
}

static PointXY rotatePlus90AboutPt(const PointXY& xy, const PointXY& origin) {
    return PointXY{-xy.y() + origin.x() + origin.y(), xy.x() - origin.x() + origin.y()};
}

static PointXY rotateMinus90AboutPt(const PointXY& xy, const PointXY& origin) {
    return PointXY{xy.y() + origin.x() - origin.y(), -xy.x() + origin.x() + origin.y()};
}

static PointXY rotatePlus180AboutPt(const PointXY& xy, const PointXY& origin) {
    return PointXY{2.0 * origin.x() - xy.x(), 2.0 * origin.y() - xy.y()};
}

}  // anonymous namespace

// constructor
LFRicCubedSphereTiles::LFRicCubedSphereTiles(const eckit::Parametrisation&) {}

std::array<std::array<double, 6>, 2> LFRicCubedSphereTiles::xy2abOffsets() const {
    return {{{0., 1., 2., 3., 0., 0.}, {1., 1., 1., 1., 2., 0.}}};
}

std::array<std::array<double, 6>, 2> LFRicCubedSphereTiles::ab2xyOffsets() const {
    return {{{0., 90., 180., 270., 0., 0.}, {-45., -45., -45., -45., 45., -135.}}};
}

void tile0Rotate(double xyz[]) {
    // Face 0, no rotation.
}

void tile1Rotate(double xyz[]) {
    double xyz_in[3];
    std::copy(xyz, xyz + 3, xyz_in);
    xyz[XX] = -xyz_in[YY];
    xyz[YY] = xyz_in[XX];
}

void tile2Rotate(double xyz[]) {
    double xyz_in[3];
    std::copy(xyz, xyz + 3, xyz_in);
    xyz[XX] = -xyz_in[XX];
    xyz[YY] = -xyz_in[YY];
}

void tile3Rotate(double xyz[]) {
    double xyz_in[3];
    std::copy(xyz, xyz + 3, xyz_in);
    xyz[XX] = xyz_in[YY];
    xyz[YY] = -xyz_in[XX];
}

void tile4Rotate(double xyz[]) {
    double xyz_in[3];
    std::copy(xyz, xyz + 3, xyz_in);
    xyz[XX] = xyz_in[ZZ];
    xyz[ZZ] = -xyz_in[XX];
}

void tile5Rotate(double xyz[]) {
    double xyz_in[3];
    std::copy(xyz, xyz + 3, xyz_in);
    xyz[XX] = -xyz_in[ZZ];
    xyz[ZZ] = xyz_in[XX];
}

void tile0RotateInverse(double xyz[]) {
    // Face 0, no rotation.
}

void tile1RotateInverse(double xyz[]) {
    double xyz_in[3];
    std::copy(xyz, xyz + 3, xyz_in);
    xyz[XX] = xyz_in[YY];
    xyz[YY] = -xyz_in[XX];
}

void tile2RotateInverse(double xyz[]) {
    double xyz_in[3];
    std::copy(xyz, xyz + 3, xyz_in);
    xyz[XX] = -xyz_in[XX];
    xyz[YY] = -xyz_in[YY];
}

void tile3RotateInverse(double xyz[]) {
    double xyz_in[3];
    std::copy(xyz, xyz + 3, xyz_in);
    xyz[XX] = -xyz_in[YY];
    xyz[YY] = xyz_in[XX];
}

void tile4RotateInverse(double xyz[]) {
    double xyz_in[3];
    std::copy(xyz, xyz + 3, xyz_in);
    xyz[XX] = -xyz_in[ZZ];
    xyz[ZZ] = xyz_in[XX];
}

void tile5RotateInverse(double xyz[]) {
    double xyz_in[3];
    std::copy(xyz, xyz + 3, xyz_in);
    xyz[XX] = xyz_in[ZZ];
    xyz[ZZ] = -xyz_in[XX];
}

void LFRicCubedSphereTiles::rotate(idx_t t, double xyz[]) const {
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
            Log::info() << "ERROR: t out of range" << std::endl;
            throw_OutOfRange("t", t, 6);
        }
    }
}

void LFRicCubedSphereTiles::unrotate(idx_t t, double xyz[]) const {
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

idx_t LFRicCubedSphereTiles::indexFromXY(const double xy[]) const {
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

    idx_t t{-1};  // tile index

    if ((xy[XX] >= 0.) && (xy[YY] >= -45.) && (xy[XX] < 90.) && (xy[YY] < 45.)) {
        t = 0;
    }
    else if ((xy[XX] >= 90.) && (xy[YY] >= -45.) && (xy[XX] < 180.) && (xy[YY] < 45.)) {
        t = 1;
    }
    else if ((xy[XX] >= 180.) && (xy[YY] >= -45.) && (xy[XX] < 270.) && (xy[YY] < 45.)) {
        t = 2;
    }
    else if ((xy[XX] >= 270.) && (xy[YY] >= -45.) && (xy[XX] < 360.) && (xy[YY] < 45.)) {
        t = 3;
    }
    else if ((xy[XX] >= 0.) && (xy[YY] >= 45.) && (xy[XX] <= 90.) && (xy[YY] <= 135.)) {
        t = 4;
    }
    else if ((xy[XX] > 0.) && (xy[YY] > -135.) && (xy[XX] < 90.) && (xy[YY] < -45.)) {
        t = 5;
    }

    return t;
}

idx_t LFRicCubedSphereTiles::indexFromLonLat(const double crd[]) const {
    idx_t t{-1};  // tile index

    double xyz[3];

    sphericalToCartesian(crd, xyz);

    const double& lon = crd[LON];

    double zPlusAbsX  = xyz[ZZ] + std::abs(xyz[XX]);
    double zPlusAbsY  = xyz[ZZ] + std::abs(xyz[YY]);
    double zMinusAbsX = xyz[ZZ] - std::abs(xyz[XX]);
    double zMinusAbsY = xyz[ZZ] - std::abs(xyz[YY]);

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

    if ((zPlusAbsX <= 0.) && (zPlusAbsY <= 0.)) {
        t = 4;
    }
    else if ((zMinusAbsX > 0.) && (zMinusAbsY > 0.)) {
        t = 5;
    }
    else if (is_geq(lon, 315.) || is_less(lon, 45.)) {
        t = 0;
    }
    else if (is_geq(lon, 45.) && is_less(lon, 135.)) {
        t = 1;
    }
    else if (is_geq(lon, 135.) && is_less(lon, 225.)) {
        t = 2;
    }
    else if (is_geq(lon, 225.) && is_less(lon, 315.)) {
        t = 3;
    }

    return t;
}

void LFRicCubedSphereTiles::enforceXYdomain(double xy[]) const {
    // the conversion from lonlat to xy can involve some small errors and small errors
    // can affect whether xy that is created within a valid space
    // This has been tested with N = 512 with equiangular and equidistant projections.
    const double tol{70.0};
    constexpr double epsilon = std::numeric_limits<double>::epsilon();

    if (debug) {
        Log::info() << "enforceXYDomain before " << xy[XX] << " " << xy[YY] << std::endl;
    }

    if (is_same(xy[YY], 45.0, tol)) {
        xy[YY] = 45.0;
    }
    if (is_same(xy[YY], -45.0, tol)) {
        xy[YY] = -45.0;
    }
    if (xy[YY] < -45.0) {
        xy[XX] = std::min(xy[XX], 90.0 - epsilon);
        xy[XX] = std::max(xy[XX], 0.0 + epsilon);
    }
    else if (xy[YY] >= -45.0) {
        xy[XX] = std::max(xy[XX], 0.0);
    }
    if (xy[YY] >= 45.0) {
        xy[XX] = std::min(xy[XX], 90.0);
    }

    xy[XX] = std::max(xy[XX], 0.0);
    xy[XX] = std::min(xy[XX], 360.0 - epsilon);
    xy[YY] = std::max(xy[YY], -135.0 + epsilon);
    xy[YY] = std::min(xy[YY], 135.0);


    if (debug) {
        Log::info() << "enforceXYDomain after " << xy[XX] << " " << xy[YY] << std::endl;
    }
}

// input is the xy value as PointXY that is a continous "cross " extension in terms of xy space from the tile in question
// the output is an xy value that lives on the standard "|---" shape
PointXY LFRicCubedSphereTiles::tileCubePeriodicity(const PointXY& xyExtended, const idx_t t) const {
    // xy space is a function of tiles--Tile 0)                    // xy space for Tile 1
    //                                                             //
    //   y ^                                                       //   y ^
    //                                                             //     |
    //     |
    //    225   ----------                                         //    225              ----------
    //     |   |To 2      |                                        //     |              |To 3 as   |
    //     |   |rotate    |                                        //     |              |below inst|
    //     |   |-180      |                                        //     |              |and rotate|
    //     |   |about     |                                        //     |              |+90 about |
    //     |   |(90,90)   |                                        //     |              |(0,45) +  |replace in [0,360) x range
    //    135  4----------4                                        //    135             4----------4
    //     |   |    <=    |                                        //     |              |To 4      |
    //     |   |          |                                        //     |              |rotate    |
    //     |   |=<   4  <=|                                        //     |              |+90 about |
    //     |   |     v    |                                        //     |              |(90,45)   |
    //     |   |*   <=    |                                        //     |              |          |
    //     45  4----------4----------4----------4----------        //     45  4----------4----------4----------4----------
    //     |   |    ^     |     ^    |    ^     |     ^    |       //     |   |    ^     |     ^    |    ^     |     ^    |
    //     |   |          |          |          |          |       //     |   |          |          |          |          |
    //     |   |=<  0    <|=<   1   <|=<  2    <|=<   3   <|       //     |   |=<  0    <|=<   1   <|=<  2    <|=<   3   <|
    //     |   |    v     |     v    |    v     |     v    |       //     |   |    v     |     v    |    v     |     v    |
    //     |   |*   =     |*    =    |*    =    | *   =    |       //     |   |*   =     |*    =    |*    =    | *   =    |
    //    -45  0 ---------1----------2----------3----------        //    -45  0 ---------1----------2----------3----------
    //     |   |    ^     |                                        //     |              |To 5      |
    //     |   |          |                                        //     |              |rotate    |
    //     |   |<   5    <|                                        //     |              |-90 about |
    //     |   |          |                                        //     |              |(90,-45)  |
    //     |   |*    v    |                                        //     |              |          |
    //   -135   ---------- ---------- ---------- ----------        //   -135   ---------- ---------- ---------- ----------
    //     ----0---------90--------180--------270--------360--->  x   ----0---------90--------180--------270--------360--->  x

    // xy space for Tile 2                                         // xy space for Tile 3
    //                                                             //
    //   y ^                                                       //   y ^
    //     |                                                       //     |
    //    225                         ----------                   //    225                                    ----------
    //     |                         |To 0      |                  //     |                                    |To 1      |
    //     |                         |rotate-180|                  //     |                                    |rotate-180|
    //     |                         |about     |                  //     |                                    |about     |
    //     |                         |(135.,90) |                  //     |                                    |(225,90)  |
    //     |                         |          |                  //     |                                    |          |
    //    135                        4----------4                  //    135                                   4----------4
    //     |                         |To 4 dble |                  //     |                                    |To 4      |
    //     |                         |+90       |                  //     |                                    |-90       |
    //     |                         |rotations |                  //     |                                    |rotation  |
    //     |                         |about(180,45)                //     |                                    |about(360,45)
    //     |                         |(90,45)   |                  //     |                                    |and       | replace in x range
    //     45  4----------4----------4----------4----------        //     45  4----------4----------4----------4----------
    //     |   |    ^     |     ^    |    ^     |     ^    |       //     |   |    ^     |     ^    |    ^     |     ^    |
    //     |   |          |          |          |          |       //     |   |          |          |          |          |
    //     |   |=<  0    <|=<   1   <|=<  2    <|=<   3   <|       //     |   |=<  0    <|=<   1   <|=<  2    <|=<   3   <|
    //     |   |    v     |     v    |    v     |     v    |       //     |   |    v     |     v    |    v     |     v    |
    //     |   |*   =     |*    =    |*    =    | *   =    |       //     |   |*   =     |*    =    |*    =    | *   =    |
    //    -45  0 ---------1----------2----------3----------        //    -45  0 ---------1----------2----------3----------
    //     |                         |To 5 dble |                  //     |                                    |To 5      |
    //     |                         |-90       |                  //     |                                    |90        |
    //     |                         |rotations |                  //     |                                    |rotation  |
    //     |                         |about (180,-45)|             //     |                                    |about (360,-45)|
    //     |                         |(90,-45)  |                  //     |                                    |and       |replace in x range
    //   -135   ---------- ---------- ---------- ----------        //   -135   ---------- ---------- ---------- ----------
    //     ----0---------90--------180--------270--------360--->  x//     ----0---------90--------180--------270--------360--->  x

    // xy space for Tile 4's                                       // xy space Tile 5's
    //                                                             //
    //   y ^                                                       //   y ^
    //     |                                                       //     |
    //    225   ----------                                         //    225   ----------
    //     |   |To 2      |                                        //     |   |To 2      |
    //     |   |rotate    |                                        //     |   |rotate    |
    //     |   |-180      |                                        //     |   |-180      |
    //     |   |about     |                                        //     |   |about     |
    //     |   |(135.,90) |                                        //     |   |(135.,90) |
    //    135  4----------4----------4----------4----------|       //    135  4----------4
    //     |   |    <=    |To 1      |To 5      |To 3      |       //     |   |    <=    |
    //     |   |          |rotate    |rotate 180|rotate+90 |       //     |   |          |
    //     |   |=<   4  <=|-90 about |about     |about     |       //     |   |=<   4  <=|
    //     |   |     v    |(90,45)   |          |          |       //     |   |     v    |
    //     |   |*   <=    |          |(135, 0)  |(360, 45) |       //     |   |*   <=    |
    //     45  4----------4----------4----------4----------        //     45  4----------4
    //     |   |    ^     |                                        //     |   |    ^     |
    //     |   |          |                                        //     |   |          |
    //     |   |=<  0    <|                                        //     |   |=   0    >|
    //     |   |    v     |                                        //     |   |          |
    //     |   |*   =     |                                        //     |   |*   =     |
    //    -45  0 ---------1                                        //    -45  0 ---------1----------2----------3----------|
    //     |   |    ^     |                                        //     |   |    ^     |To 1      |To 4      |To 3 as   |
    //     |   |          |                                        //     |   |          |rotate    |rotate    |left instr|
    //     |   |<   5    <|                                        //     |   |<   5    <|+90 about |+180      |+ rotate  |
    //     |   |          |                                        //     |   |          |(90,-45)  |(135,0)   |+90 about |
    //     |   |*    v    |                                        //     |   |*    v    |          |          |(0, 45)   |+ replace in [0,360 range)
    //   -135   ---------- ---------- ---------- ----------        //   -135   ---------- ---------- ---------- ----------
    //     ----0---------90--------180--------270--------360--->  x//     ----0---------90--------180--------270--------360--->  x


    // Step 1: wrap into range  x = [ 0, 360],  y = [135, 225] ensuring that we stay
    //         on the cross associated with tile index.
    if (!withinCross(t, xyExtended)) {
        ATLAS_THROW_EXCEPTION("Point (" << xyExtended.x() << "," << xyExtended.y() << ")"
                                        << " is not in the cross extension of tile " << t);
    }

    PointXY withinRange = xyExtended;
    enforceWrapAround(t, withinRange);

    PointXY finalXY = withinRange;
    PointXY tempXY  = withinRange;

    switch (t) {
        case 0:
            finalXY = (withinRange.y() > 135.0) ? rotatePlus180AboutPt(withinRange, PointXY{135.0, 90.0}) : withinRange;
            break;
        case 1:
            if ((withinRange.x() >= 90.0) && (withinRange.x() <= 180.0)) {
                if (withinRange.y() >= 45.0) {
                    tempXY = rotatePlus90AboutPt(withinRange, PointXY{90.0, 45.0});

                    if (withinRange.y() > 135.0) {
                        finalXY = rotatePlus90AboutPt(tempXY, PointXY{0.0, 45.0});
                        finalXY.x() += 360.0;
                    }
                    else {
                        finalXY = tempXY;
                    }
                }
                else if (withinRange.y() < -45.0) {
                    finalXY = rotateMinus90AboutPt(withinRange, PointXY{90.0, -45.0});
                }
            }
            break;
        case 2:
            if ((withinRange.x() >= 180.0) && (withinRange.x() <= 270.0)) {
                if (withinRange.y() > 135.0) {
                    finalXY = rotatePlus180AboutPt(tempXY, PointXY{135.0, 90.0});
                }
                else if (withinRange.y() >= 45.0) {
                    tempXY  = rotatePlus90AboutPt(withinRange, PointXY{180.0, 45.0});
                    finalXY = rotatePlus90AboutPt(tempXY, PointXY{90.0, 45.0});
                }
                else if (withinRange.y() < -45.0) {
                    tempXY  = rotateMinus90AboutPt(withinRange, PointXY{180.0, -45.0});
                    finalXY = rotateMinus90AboutPt(tempXY, PointXY{90.0, -45.0});
                }
            }
            break;

        case 3:
            if (((withinRange.x() >= 270.0) && (withinRange.x() <= 360.0)) || withinRange.x() == 0.0) {
                if (withinRange.y() > 135.0) {
                    finalXY = rotatePlus180AboutPt(tempXY, PointXY{225.0, 90.0});
                }
                else if (withinRange.y() >= 45.0) {
                    finalXY = rotateMinus90AboutPt(tempXY, PointXY{360.0, 45.0});
                    finalXY.x() -= 360.0;
                }
                else if (withinRange.y() < -45.0) {
                    finalXY = rotatePlus90AboutPt(tempXY, PointXY{360.0, -45.0});
                    finalXY.x() -= 360.0;
                }
            }
            break;
        case 4:
            if (withinRange.y() > 135.0) {
                finalXY = rotatePlus180AboutPt(tempXY, PointXY{135.0, 90.0});
            }
            else if ((withinRange.y() >= 45.0) && (withinRange.y() <= 135.0)) {
                if ((withinRange.x() > 90.0) && (withinRange.x() <= 180.0)) {
                    finalXY = rotateMinus90AboutPt(withinRange, PointXY{90.0, 45.0});
                }
                if ((withinRange.x() > 180.0) && (withinRange.x() <= 270.0)) {
                    finalXY = rotatePlus180AboutPt(withinRange, PointXY{135.0, 0.0});
                }
                if ((withinRange.x() > 270.0) && (withinRange.x() <= 360.0)) {
                    finalXY = rotatePlus90AboutPt(withinRange, PointXY{360.0, 45.0});
                }
            }
            break;
        case 5:
            if (withinRange.y() > 135.0) {
                finalXY = rotatePlus180AboutPt(tempXY, PointXY{135.0, 90.0});
            }
            else if ((withinRange.y() <= -45.0) && (withinRange.y() >= -135.0)) {
                if ((withinRange.x() > 90.0) && (withinRange.x() <= 180.0)) {
                    finalXY = rotatePlus90AboutPt(withinRange, PointXY{90.0, -45.0});
                }
                if ((withinRange.x() > 180.0) && (withinRange.x() <= 270.0)) {
                    finalXY = rotatePlus180AboutPt(withinRange, PointXY{135.0, 0.0});
                }
                if ((withinRange.x() > 270.0) && (withinRange.x() <= 360.0)) {
                    finalXY = rotateMinus90AboutPt(withinRange, PointXY{360.0, -45.0});
                }
            }
            break;
    }

    // Now we are on the standard tiles with their default rotate orientation.
    // The next step is to wrap the edges and corners to the correct tiles.
    if (finalXY.x() > 90. && finalXY.x() < 180. && finalXY.y() == 45.0) {
        finalXY.y() = 45.0 + (finalXY.x() - 90.0);
        finalXY.x() = 90.0;
    }
    if (finalXY.x() == 180. && finalXY.y() == 45.0) {
        finalXY.x() = 90.;
        finalXY.y() = 135.0;
    }
    if (finalXY.x() > 180. && finalXY.x() < 270. && finalXY.y() == 45.0) {
        finalXY.x() = 90.0 - (finalXY.x() - 90.0);
        finalXY.y() = 135.0;
    }
    if (finalXY.x() == 270. && finalXY.y() == 45.0) {
        finalXY.x() = 0.;
        finalXY.y() = 135.0;
    }
    if (finalXY.x() > 270. && finalXY.x() <= 360. && finalXY.y() == 45.0) {
        finalXY.y() = 135.0 - (finalXY.x() - 270.0);
        finalXY.x() = 0.0;
    }

    if (finalXY.x() >= 360.) {
        finalXY.x() -= 360.;
    }

    if (finalXY.x() == 0. && finalXY.y() < -45. && finalXY.y() > -135.0) {
        finalXY.x() = 360.0 + (finalXY.y() + 45.0);
        finalXY.y() = -45.;
    }
    if (finalXY.x() == 90. && finalXY.y() < -45. && finalXY.y() > -135.0) {
        finalXY.x() = 90.0 - (finalXY.y() + 45.0);
        finalXY.y() = -45.;
    }
    if (finalXY.y() <= -135. && finalXY.x() >= 0. && finalXY.x() <= 90.0) {
        finalXY.x() = 270.0 - finalXY.x();
        finalXY.y() = -45.;
    }

    return finalXY;
}


void LFRicCubedSphereTiles::print(std::ostream& os) const {
    os << "CubedSphereTiles["
       << "]";
}


bool LFRicCubedSphereTiles::withinCross(const idx_t tiles, const PointXY& withinRange) const {
    std::size_t t = static_cast<std::size_t>(tiles);
    return !((withinRange.x() < botLeftTile(t).x() && withinRange.y() < botLeftTile(t).y()) ||
             (withinRange.x() > botRightTile(t).x() && withinRange.y() < botRightTile(t).y()) ||
             (withinRange.x() > topRightTile(t).x() && withinRange.y() > topRightTile(t).y()) ||
             (withinRange.x() < topLeftTile(t).x() && withinRange.y() > topLeftTile(t).y()));
}


void LFRicCubedSphereTiles::enforceWrapAround(const idx_t t, PointXY& withinRange) const {
    if (withinRange.x() < 0.0) {
        PointXY temp = withinRange;
        temp.x() += 360;
        if (withinCross(t, temp)) {
            withinRange = temp;
        }
    }
    if (withinRange.x() > 360.0) {
        PointXY temp = withinRange;
        temp.x() -= 360;
        if (withinCross(t, temp)) {
            withinRange = temp;
        }
    }
    if (withinRange.y() <= -135.0) {
        PointXY temp = withinRange;
        temp.y() += 360;
        if (withinCross(t, temp)) {
            withinRange = temp;
        }
    }
    if (withinRange.y() > 225.0) {
        PointXY temp = withinRange;
        temp.y() -= 360;
        if (withinCross(t, temp)) {
            withinRange = temp;
        }
    }

    return;
}

const PointXY& LFRicCubedSphereTiles::tileCentre(size_t t) const {
    return tileCentres_[t];
}

const LFRicCubedSphereTiles::Jacobian& LFRicCubedSphereTiles::tileJacobian(size_t t) const {
    return tileJacobians_[t];
}

PointXY LFRicCubedSphereTiles::botLeftTile(size_t t) {
    return tileCentres_[t] + PointXY{-45., -45.};
}

PointXY LFRicCubedSphereTiles::botRightTile(size_t t) {
    return tileCentres_[t] + PointXY{45., -45.};
}

PointXY LFRicCubedSphereTiles::topLeftTile(size_t t) {
    return tileCentres_[t] + PointXY{-45., 45.};
}

PointXY LFRicCubedSphereTiles::topRightTile(size_t t) {
    return tileCentres_[t] + PointXY{45., 45.};
}


// Centre of each tile in xy space.
const std::array<PointXY, 6> LFRicCubedSphereTiles::tileCentres_{
    PointXY{45., 0.}, PointXY{135., 0.}, PointXY{225., 0.}, PointXY{315., 0.}, PointXY{45., 90.}, PointXY{45., -90.}};

// Jacobian of xy space with respect to curvilinear coordinates for each tile.
const std::array<LFRicCubedSphereTiles::Jacobian, 6> LFRicCubedSphereTiles::tileJacobians_{
    Jacobian{{1., 0.}, {0., 1.}},  Jacobian{{1., 0.}, {0., 1.}}, Jacobian{{0., -1.}, {1., 0.}},
    Jacobian{{0., -1.}, {1., 0.}}, Jacobian{{1., 0.}, {0., 1.}}, Jacobian{{0., 1.}, {-1., 0.}}};


namespace {
static CubedSphereTilesBuilder<LFRicCubedSphereTiles> register_builder(LFRicCubedSphereTiles::static_type());
}


}  // namespace detail
}  // namespace grid
}  // namespace atlas
