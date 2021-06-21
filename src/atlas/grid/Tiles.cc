/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>

#include "atlas/grid/Tiles.h"
#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/grid/detail/tiles/FV3Tiles.h"
#include "atlas/grid/detail/tiles/LFRicTiles.h"


using FV3CubedSphereTiles = atlas::cubedspheretiles::FV3CubedSphereTiles;
using LFRicCubedSphereTiles = atlas::cubedspheretiles::LFRicCubedSphereTiles;

namespace atlas {

CubedSphereTiles::CubedSphereTiles( const eckit::Parametrisation& p ) : Handle(
                atlas::cubedspheretiles::CubedSphereTiles::create( p ) ) {

   std::cout << "grid/Tiles.h CubedSphereTiles constr"  << std::endl;
}

FV3CubedSphereTiles::FV3CubedSphereTiles( const CubedSphereTiles& cubedspheretiles ) :
    CubedSphereTiles ( cubedspheretiles ),
    cubedspheretiles_( dynamic_cast<const atlas::cubedspheretiles::FV3CubedSphereTiles*>( get() ) ) {}


LFRicCubedSphereTiles::LFRicCubedSphereTiles( const CubedSphereTiles& cubedspheretiles ) :
    CubedSphereTiles ( cubedspheretiles ),
    cubedspheretiles_( dynamic_cast<const atlas::cubedspheretiles::LFRicCubedSphereTiles*>( get() ) ) {}


std::string atlas::CubedSphereTiles::type() const {
    return get()->type();
}

void CubedSphereTiles::tile0Rotate( double xyz[] ) const {
    return get()->tile0Rotate(xyz);
}

void CubedSphereTiles::tile1Rotate( double xyz[] ) const{
    return get()->tile1Rotate(xyz);
};

void CubedSphereTiles::tile2Rotate( double xyz[] ) const{
    return get()->tile2Rotate(xyz);
};

void CubedSphereTiles::tile3Rotate( double xyz[] ) const {
    return get()->tile3Rotate(xyz);
}

void CubedSphereTiles::tile4Rotate( double xyz[] ) const {
    return get()->tile4Rotate(xyz);
}

void CubedSphereTiles::tile5Rotate( double xyz[] ) const{
    return get()->tile5Rotate(xyz);
}

void CubedSphereTiles::tile0RotateInverse( double xyz[] ) const {
    return get()->tile0RotateInverse(xyz);
}

void CubedSphereTiles::tile1RotateInverse( double xyz[] ) const {
    return get()->tile1RotateInverse(xyz);
}

void CubedSphereTiles::tile2RotateInverse( double xyz[] ) const {
    return get()->tile2RotateInverse(xyz);
}

void CubedSphereTiles::tile3RotateInverse( double xyz[] ) const {
    return get()->tile3RotateInverse(xyz);
}

void CubedSphereTiles::tile4RotateInverse( double xyz[] ) const {
    return get()->tile4RotateInverse(xyz);
}

void CubedSphereTiles::tile5RotateInverse( double xyz[] ) const {
    return get()->tile5RotateInverse(xyz);
}

idx_t CubedSphereTiles::tileFromXY( const double xy[] ) const {
     return get()->tileFromXY(xy);
}

idx_t CubedSphereTiles::tileFromLonLat(const double lonlat[]) const {
     return get()->tileFromLonLat(lonlat);
}

void CubedSphereTiles::enforceXYdomain(double xy[]) const {
     return get()->enforceXYdomain(xy);
}

void CubedSphereTiles::print( std::ostream& os ) const {
    get()->print( os );
}

std::ostream& operator<<( std::ostream& os, const CubedSphereTiles& t ) {
    t.print( os );
    return os;
}

}  // namespace atlas
