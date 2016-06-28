/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/global/lonlat/LonLat.h"


namespace atlas {
namespace grid {
namespace global {
namespace lonlat {


std::string LonLat::grid_type_str() {
    return "reduced_lonlat";
}


std::string LonLat::className() {
    return "atlas.grid.global.lonlat.LonLat";
}


LonLat::LonLat(const Shift& shift, const Domain& dom) :
    Structured(dom),
    shift_(shift) {
}


} // namespace lonlat
} // namespace global
} // namespace grid
} // namespace atlas
