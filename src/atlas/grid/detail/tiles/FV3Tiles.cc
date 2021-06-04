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
#include <string>
#include <utility>

#include "eckit/utils/Hash.h"

#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/grid/detail/tiles/FV3Tiles.h"

namespace atlas {
namespace cubedspheretiles {

// constructor
FV3CubedSphereTiles::FV3CubedSphereTiles( const eckit::Parametrisation& ) {

}

idx_t FV3CubedSphereTiles::tileFromXY( const double xy[] ) const  {
    idx_t t{-1}; // tile index

    return t;
}

idx_t FV3CubedSphereTiles::tileFromLonLat( const double lonlat[] ) const {
    idx_t t{-1}; // tile index

    return t;
}

void FV3CubedSphereTiles::enforceXYdomain( double xy[] ) const {

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
