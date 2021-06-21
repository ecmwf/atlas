/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string>
#include <iostream>

#include "eckit/utils/MD5.h"

#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/grid/detail/tiles/TilesFactory.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace cubedspheretiles {

const CubedSphereTiles* CubedSphereTiles::create() {
    // default: FV3 version (for now)
    util::Config params;
    params.set( "tile type", "FV3CubedSphereTiles" );
    return CubedSphereTiles::create( params );
}

const CubedSphereTiles* CubedSphereTiles::create( const eckit::Parametrisation& p ) {
    std::string CubedSphereTiles_type;

    std::cout << "tiles/Tiles.cc within create " << std::endl;
    if ( p.get( "tile type", CubedSphereTiles_type ) ) {
        std::cout << "tiles/Tiles.cc tile type = "  << CubedSphereTiles_type << std::endl;
        return CubedSphereTilesFactory::build( CubedSphereTiles_type, p );
    } else {
        return create();
    }

}




}  // namespace cubedspheretiles
}  // namespace atlas
