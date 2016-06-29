/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <typeinfo>
#include "eckit/memory/Builder.h"
#include "atlas/internals/atlas_config.h"
#include "atlas/grid/global/lonlat/ShiftedLat.h"


namespace atlas {
namespace grid {
namespace global {
namespace lonlat {


register_BuilderT1(Grid,ShiftedLat,ShiftedLat::grid_type_str());


std::string ShiftedLat::grid_type_str() {
    return "shifted_lat";
}


std::string ShiftedLat::className() {
    return "atlas.grid.global.lonlat.ShiftedLat";
}


void ShiftedLat::set_typeinfo() {
    std::stringstream s;
    if( N() ) {
        s << "Slat" << N();
    } else {
        s << "Slat" << nlon() << "x" << nlat();
    }
    shortName_ = s.str();
    grid_type_ = grid_type_str();
}


ShiftedLat::ShiftedLat(const eckit::Parametrisation& p) :
    LonLat(Shift::LAT,Domain::makeGlobal()) {
    setup(p);
    set_typeinfo();
}


ShiftedLat::ShiftedLat(const size_t N) :
    LonLat(Shift::LAT,Domain::makeGlobal()) {
    LonLat::setup(N,domain());
}


ShiftedLat::ShiftedLat(const size_t nlon, const size_t nlat) :
    LonLat(Shift::LAT,Domain::makeGlobal()) {
    LonLat::setup(nlon,nlat,domain());
    set_typeinfo();
}


void ShiftedLat::setup(const eckit::Parametrisation& p) {
    // set nlon and nlat
    size_t nlon, nlat, N(0);

    if (p.get("N",N)) {
        LonLat::setup(N,domain_);
    }
    else if (p.get("nlon", nlon) && p.get("nlat", nlat)) {
        LonLat::setup(nlon, nlat, domain_);
    }
    else {
        throw eckit::BadParameter("Params (nlon,nlat) or N missing", Here());
    }
}


eckit::Properties ShiftedLat::spec() const {
    eckit::Properties grid_spec;

    grid_spec.set("grid_type",gridType() );
    grid_spec.set("short_name",shortName());

    grid_spec.set("N", N() );
    grid_spec.set("nlon", nlon() );
    grid_spec.set("nlat", nlat() );

    return grid_spec;
}


extern "C" {


    Structured* atlas__grid__global__lonlat__ShiftedLat(size_t nlon, size_t nlat) {
        return new ShiftedLat(nlon,nlat);
    }


}


} // namespace global
} // namespace lonlat
} // namespace grid
} // namespace atlas
