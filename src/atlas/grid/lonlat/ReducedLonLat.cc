/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/lonlat/ReducedLonLat.h"


namespace atlas {
namespace grid {
namespace lonlat {


register_BuilderT1(Grid,ReducedLonLat,ReducedLonLat::grid_type_str());


std::string ReducedLonLat::grid_type_str() {
    return "reduced_lonlat";
}


std::string ReducedLonLat::className() {
    return "atlas.grid.lonlat.ReducedLonLat";
}


void ReducedLonLat::set_typeinfo() {
    std::stringstream s;
    s << "reduced_lonlat";
    shortName_ = s.str();
    grid_type_ = grid_type_str();
}


ReducedLonLat::ReducedLonLat(const size_t nlat, const long nlon[], const Domain& dom) :
    LonLat(Shift::NONE,dom) {
    LonLat::setup(nlon,nlat,domain());
    set_typeinfo();
}


ReducedLonLat::ReducedLonLat(const eckit::Parametrisation& params) :
    LonLat(Shift::NONE,Domain::makeGlobal()) {
    setup(params);
    set_typeinfo();
}


void ReducedLonLat::setup(const eckit::Parametrisation& params) {
    if (!params.has("nlat"))  throw eckit::BadParameter("N missing in Params",Here());
    if (!params.has("pl"))    throw eckit::BadParameter("npts_per_lat missing in Params",Here());

    size_t nlat;
    std::vector<long> nlon;
    params.get("nlat", nlat);
    params.get("pl", nlon);
    params.get("N", Structured::N_);

    LonLat::setup(nlon.data(), nlat, domain());

    ASSERT(shift_(Shift::NONE));
}


eckit::Properties ReducedLonLat::spec() const {
    eckit::Properties grid_spec;

    grid_spec.set("grid_type",  gridType());
    grid_spec.set("short_name", shortName());
    grid_spec.set("N",          N());
    grid_spec.set("nlat",       nlat());
    grid_spec.set("pl",         eckit::makeVectorValue(pl()));
    grid_spec.set("latitudes",  eckit::makeVectorValue(latitudes()));

    return grid_spec;
}


}  // namespace lonlat
}  // namespace grid
}  // namespace atlas

