/*
 * (C) Crown Copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string>

#include "atlas/grid/detail/tiles/FV3Tiles.h"
#include "atlas/grid/detail/tiles/LFRicTiles.h"
#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/grid/detail/tiles/TilesFactory.h"
#include "atlas/util/Factory.h"

namespace atlas {
namespace grid {
namespace detail {

//----------------------------------------------------------------------------------------------------------------------

namespace {
void force_link() {
    static struct Link {
        Link() {
            CubedSphereTilesBuilder<FV3CubedSphereTiles>();
            CubedSphereTilesBuilder<LFRicCubedSphereTiles>();
        }
    } link;
}
}  // namespace

//----------------------------------------------------------------------------------------------------------------------

const CubedSphereTiles* CubedSphereTilesFactory::build(const std::string& builder) {
    return build(builder, util::NoConfig());
}

const CubedSphereTiles* CubedSphereTilesFactory::build(const std::string& builder,
                                                       const eckit::Parametrisation& param) {
    force_link();
    auto factory = get(builder);
    return factory->make(param);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace grid
}  // namespace atlas
