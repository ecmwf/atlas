/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/grid/detail/tiles/LFRicTiles.h"

namespace atlas {
namespace cubedspheretiles {


// constructor
LFRicCubedSphereTiles::LFRicCubedSphereTiles( const eckit::Parametrisation& ) {

}

idx_t LFRicCubedSphereTiles::tileFromXY( const double xy[] ) const  {

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

    idx_t t{-1}; // tile index

    return t;
}

idx_t LFRicCubedSphereTiles::tileFromLonLat( const double lonlat[] ) const {
    idx_t t{-1}; // tile index

    return t;
}

void LFRicCubedSphereTiles::enforceXYdomain( double xy[] ) const {

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
