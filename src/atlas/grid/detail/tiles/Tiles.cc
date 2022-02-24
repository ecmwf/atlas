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
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {
namespace detail {

const CubedSphereTiles* CubedSphereTiles::create() {
    // default: LFRic version (for now)
    util::Config params;
    params.set("type", "cubedsphere_lfric");
    return CubedSphereTiles::create(params);
}

const CubedSphereTiles* CubedSphereTiles::create(const std::string& s) {
    util::Config params;
    if (s == "") {
        return CubedSphereTiles::create();
    }
    params.set("type", s);
    return CubedSphereTiles::create(params);
}

const CubedSphereTiles* CubedSphereTiles::create(const eckit::Parametrisation& p) {
    std::string CubedSphereTiles_type;

    if (p.has("type")) {
        p.get("type", CubedSphereTiles_type);
        return CubedSphereTilesFactory::build(CubedSphereTiles_type, p);
    }
    else {
        return create();
    }
}

}  // namespace detail
}  // namespace grid
}  // namespace atlas
