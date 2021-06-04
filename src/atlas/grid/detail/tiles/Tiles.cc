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

#include "eckit/utils/MD5.h"

#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/grid/detail/tiles/TilesFactory.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace cubedspheretiles {

const CubedSphereTiles* CubedSphereTiles::create( const eckit::Parametrisation& p ) {
    std::string CubedSphereTiles_type;
    if ( p.get( "type", CubedSphereTiles_type ) ) {
        return CubedSphereTilesFactory::build( CubedSphereTiles_type, p );
    }

    // should return error here
    throw_Exception( "type missing in Params", Here() );
}

}  // namespace cubedspheretiles
}  // namespace atlas
