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
