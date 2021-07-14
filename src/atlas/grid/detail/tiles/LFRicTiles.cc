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
#include <ostream>
#include <iostream>

#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/grid/detail/tiles/TilesFactory.h"
#include "atlas/grid/detail/tiles/LFRicTiles.h"
#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Log.h"
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

enum dirn {
    E,
    S,
    W,
    N
};

}

// constructor
LFRicCubedSphereTiles::LFRicCubedSphereTiles( const eckit::Parametrisation& ) {
}


// enum E S W N
/*
define node

static std::array<std::array<std::pair<int, int>, 4>, 6> tileEdgeMapping() {
    return {
        { { (3, )   }


        }
    };
}
[t][e]pair(t2, e)

[t][c]  -> t c  t c
*/

std::array<std::array<double,6>,2> LFRicCubedSphereTiles::xy2abOffsets() const {
    return { {  {0., 1., 2., 3., 0., 0.},
                {1., 1., 1., 1., 2., 0.} } };
}

std::array<std::array<double,6>,2> LFRicCubedSphereTiles::ab2xyOffsets() const {
    return { { {0.,  90., 180., 270.,  0.,    0.},
             {-45., -45., -45., -45., 45., -135.} } };
}

void LFRicCubedSphereTiles::tile0Rotate( double xyz[] ) const {
    // Face 0, no rotation.
}

void LFRicCubedSphereTiles::tile1Rotate( double xyz[] ) const {
    double xyz_in[3];
    std::copy( xyz, xyz + 3, xyz_in );
    xyz[XX] = -xyz_in[YY];
    xyz[YY] =  xyz_in[XX];
}

void LFRicCubedSphereTiles::tile2Rotate( double xyz[] ) const {
    double xyz_in[3];
    std::copy( xyz, xyz + 3, xyz_in );
    xyz[XX] = -xyz_in[XX];
    xyz[YY] = -xyz_in[YY];

}
void LFRicCubedSphereTiles::tile3Rotate( double xyz[] ) const {
    double xyz_in[3];
    std::copy( xyz, xyz + 3, xyz_in );
    xyz[XX] =  xyz_in[YY];
    xyz[YY] = -xyz_in[XX];
}
void LFRicCubedSphereTiles::tile4Rotate( double xyz[] ) const {
    double xyz_in[3];
    std::copy( xyz, xyz + 3, xyz_in );
    xyz[XX] =  xyz_in[ZZ];
    xyz[ZZ] = -xyz_in[XX];
}
void LFRicCubedSphereTiles::tile5Rotate( double xyz[] ) const {
    double xyz_in[3];
    std::copy( xyz, xyz + 3, xyz_in );
    xyz[XX] = -xyz_in[YY];
    xyz[YY] =  xyz_in[ZZ];
    xyz[ZZ] =  xyz_in[XX];

}

void LFRicCubedSphereTiles::tile0RotateInverse( double xyz[] ) const {
    // Face 0, no rotation.
}

void LFRicCubedSphereTiles::tile1RotateInverse( double xyz[] ) const {
    double xyz_in[3];
    std::copy( xyz, xyz + 3, xyz_in );
    xyz[XX] =  xyz_in[YY];
    xyz[YY] = -xyz_in[XX];
}

void LFRicCubedSphereTiles::tile2RotateInverse( double xyz[] ) const {
    double xyz_in[3];
    std::copy( xyz, xyz + 3, xyz_in );
    xyz[XX] = -xyz_in[XX];
    xyz[YY] = -xyz_in[YY];
}

void LFRicCubedSphereTiles::tile3RotateInverse( double xyz[] ) const {
    double xyz_in[3];
    std::copy( xyz, xyz + 3, xyz_in );
    xyz[XX] = -xyz_in[YY];
    xyz[YY] =  xyz_in[XX];
}

void LFRicCubedSphereTiles::tile4RotateInverse( double xyz[] ) const {
    double xyz_in[3];
    std::copy( xyz, xyz + 3, xyz_in );
    xyz[XX] = -xyz_in[ZZ];
    xyz[ZZ] =  xyz_in[XX];
}

void LFRicCubedSphereTiles::tile5RotateInverse( double xyz[] ) const {
    double xyz_in[3];
    std::copy( xyz, xyz + 3, xyz_in );
    xyz[XX] =  xyz_in[ZZ];
    xyz[YY] = -xyz_in[XX];
    xyz[ZZ] =  xyz_in[YY];
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
    } else if (lon >= 315.  || lon < 45.) {
         t = 0;
    } else if (lon >= 45. && lon < 135.) {
         t = 1;
    } else if (lon >= 135.  && lon < 225.) {
         t = 2;
    } else if (lon >= 225.  && lon < 315.) {
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
        Log::info() << "enforceXYDomain after " << xy[XX] << " " << xy[YY] << std::endl;
    }
}

// input is the xy value as PointXY that is not quite on the right part of the XY space
// i.e. a ghost point X
atlas::PointXY LFRicCubedSphereTiles::anyXYToFundamentalXY (const atlas::PointXY & anyXY ) const {


    // I think that mathematically the best option is to consider
    // this as a rotation of a tile and a translation


    // tile values are periodicity to them
    //  x = [0, 360)   y = (-135, 225]
    //  with the exception that all edges below  y=-45 get mapped onto y=-45.


    // Assume one face-edge is of length 90 degrees.
    //
    //   y ^
    //     |
    //    225  4----------4----------4----------4----------
    //     |   |    ^    *|     ^   *|    ^    *|     ^   *|
    //     |   |    =     |     =    |          |          |
    //     |   |<   2   =<|<   3   =<|<  0    =<|=<   1   <|
    //     |   |          |     v    |    v     |     v    |
    //     |   |    <     |          |     =    |     =    |
    //    135   ---------- ---------- ---------- ----------
    //     |   |    <=    |          |         *|*         |
    //     |   |          |          |          |          |
    //     |   |=<   4  <=|     4    |     4    |     4    |
    //     |   |     v    |     ^    |     ^    |     ^    |
    //     |   |*   <=    |     |   *|     |    |     |    |
    //     45  4----------4----------4----------4----------
    //     |   |    ^     |     ^    |    ^     |     ^    |
    //     |   |          |          |          |          |
    //     |   |=<  0    <|=<   1   <|=<  2    <|=<   3   <|
    //     |   |    v     |     v    |    v     |     v    |
    //     |   |*   =     |*    =    |*    =    | *   =    |
    //    -45  0 ---------1----------2----------3----------
    //     |   |    ^     |*    ^    |     ^   *|    ^     |
    //     |   |          |          |          |          |
    //     |   |<   5    <|<   5    <|<   5    <|<   5    <|
    //     |   |          |          |          |          |
    //     |   |*    v    |     v    |     v    |     v   *|
    //   -135   ---------- ---------- ---------- ----------
    //     ----0---------90--------180--------270--------360--->  x


    // first find the bottom left hand corner of xyGhost tile
    // (except for tile 5) this point will exist in the tile.
    atlas::PointXY withinRange = anyXY;

    while (withinRange.x() < 0.0) { withinRange.x() += 360.0; }
    while (withinRange.x() >= 360.0) { withinRange.x() -= 360.0; }
    while (withinRange.y() <= -135.0) { withinRange.y() += 360.0; }
    while (withinRange.y() > 225.0) { withinRange.y() -= 360.0; }

    std::cout << "withinRange " << withinRange.x() << " " << withinRange.y() << std::endl;

    // Note that the tile index here is that for the inner part of the tile
    // It won't be always consistent with the official tile indexing
    // However, it does not need to, here.
    std::size_t tileValue[4][4] = { {5, 5, 5, 5},
                                    {0, 1, 2, 3},
                                    {4, 4, 4, 4},
                                    {2, 3, 0, 1} };

    // originCorner shows the origin for the xy displacement
    // 0 = bottom left
    // 1 = bottom right
    // 2 = top right
    // 3 = top left.
    std::size_t originCorner[4][4] = { {0, 3, 2, 1 },
                                       {0, 0, 0, 0 },
                                       {0, 1, 2, 3 },
                                       {2, 2, 2, 2 } };

    std::array<atlas::PointXY, 6> botLeftTile{atlas::PointXY{0., -45.},   atlas::PointXY{90, -45},
                                              atlas::PointXY{180., -45.}, atlas::PointXY{270, -45},
                                              atlas::PointXY{0., 45.},    atlas::PointXY{0, -135.} };

    // Step 0: Find appropriate panel

    auto xIndex = static_cast<size_t>(withinRange.x()/90.0);
    auto yIndex = static_cast<size_t>((withinRange.y()+135.0)/90.0);

    // The inverted L edge of tile 4 needs including in tile 4
    if (((withinRange.x() == 90.0) || withinRange.x() == 180.0 ||  withinRange.x() == 270.0)
         && (withinRange.y() > 45.0)) xIndex -= 1;
    if (withinRange.y() == 135.0) yIndex -= 1;

    std::cout << "X Y index = " << xIndex << " " << yIndex << std::endl;

    atlas::PointXY BotLeft{ xIndex * 90.0, yIndex * 90.0 - 135.0};
    atlas::PointXY remBotLeft = withinRange - BotLeft;
    atlas::PointXY rottransPt;

    if (yIndex != 4) {

        std::size_t tile = tileValue[yIndex][xIndex];

        std::cout << "tile " << tile << std::endl;

        // Step 2: Apply appropriate rotation:
        atlas::PointXY remBotRight = atlas::PointXY{ 90.0 - remBotLeft.y(), remBotLeft.x()};
        atlas::PointXY remTopRight = atlas::PointXY{ 90.0 - remBotLeft.x(), 90.0 - remBotLeft.y()};
        atlas::PointXY remTopLeft = atlas::PointXY{ remBotLeft.y(), 90 - remBotLeft.x()};
        std::array<atlas::PointXY, 4> orientatedDisp{remBotLeft, remBotRight, remTopRight, remTopLeft};

        atlas::PointXY rem = orientatedDisp[ originCorner[yIndex][xIndex] ];

        std::cout << "rem botLeftTile " << rem << " "  << botLeftTile[tile]  << std::endl;

        rottransPt = rem + botLeftTile[tile];


    } else  {
        // yIndex equal to 4 means that you are on the top most horizontal line.
        rottransPt.x() = (withinRange.x() < 180 ?   withinRange.x() + 180.:  withinRange.x() - 180.);
        rottransPt.y() = - 45.0;

    }

    // so far we have dealt not explicitly dealt with edges
    // the only edges that need considering are those that concern the bottom tile
    // The rest, I think are either dealt with by the operation to wrap (xy) within the x = [0,360)
    // y = [-135, 215] ranges.

    if ((rottransPt.x() == 0.0) && (rottransPt.y() < -45.0)) {
       // x = [270,360)
       rottransPt.x() = 360.0 + (45.0 + rottransPt.y());
       rottransPt.y() = -45.0;
    }
    if ((rottransPt.x() == 90.0) && (rottransPt.y() < -45.0)) {
       // x = (90,180]
       rottransPt.x() = 90.0 - (45.0 + rottransPt.y());
       rottransPt.y() = -45.0;
    }
    if ((rottransPt.x() >= 0.0) && (rottransPt.x() <= 90.0) && (rottransPt.y() == -135.0)) {
       rottransPt.x() = 90.0 - rottransPt.x();
       rottransPt.y() = -45.0;
    }
    return rottransPt;
}


void LFRicCubedSphereTiles::print( std::ostream& os) const {
    os << "CubedSphereTiles["
       << "]";
}

namespace {
static  CubedSphereTilesBuilder<LFRicCubedSphereTiles> register_builder( LFRicCubedSphereTiles::static_type() );
}


/*

// I am first going to stipulate that this will work for (x = [-90, 450), y= (-225,225])
atlas::PointXY origin({0., 0.});

atlas::PointXY botLeft({ static_cast<int>(xyGhost.x()/90.0) * 90.0,
                         static_cast<int>((xyGhost.y()+135.0)/90.0) * 90.0 - 135.0
                       });
atlas::PointXY remBotLeft = xyGhost - botLeft;



// Deal with data on tiles 0-4 where botLeft is part of the tile. (L shaped edges.)
std::vector<atlas::PointXY> ownedBotLeft{ {0., -45.}, {90., -45.}, {180., -45.}, {270.0, -45.0}, {0, 45.} };
for (std::size_t i = 0; i < ownedBotLeft.size(); ++i) {
    if ((botLeft == ownedBotLeft[i]) && (remBotLeft.x() >= 0.) && (remBotLeft.y() >= 0.)) {
        return xyGhost;
    }
}

// Deal with __
//            |  shaped edges
atlas::PointXY topRight = botLeft + atlas::PointXY{90.0, 90.0};
atlas::PointXY remTopRight = topRight - xyGhost;
std::vector<atlas::PointXY> ownedTopRight{ {90., 135.}};
for (std::size_t i = 0; i < ownedTopRight.size(); ++i) {
    if ((botLeft == ownedTopRight[i]) && (remTopRight.x() >= 0.) && (remTopRight.y() >= 0.)) {
        return xyGhost;
    }
}

// Deal tiles that have no owned edges.
std::vector<atlas::PointXY> notOwnedBotLeft{ {0., -135.}};
for (std::size_t i = 0; i < notOwnedBotLeft.size(); ++i) {
    if ((botLeft == notOwnedBotLeft[i]) && (remBotLeft.x() >= 0.) && (remBotLeft.y() >= 0.) &&
        (remTopRight.x() >= 0.) && (remTopRight.y() >= 0.)) {
        return xyGhost;
    }
}

*/
}  // namespace cubespheretiles
}  // namespace atlas
