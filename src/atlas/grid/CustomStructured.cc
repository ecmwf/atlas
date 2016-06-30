/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/CustomStructured.h"


namespace atlas {
namespace grid {


register_BuilderT1(Grid, CustomStructured, CustomStructured::grid_type_str());


std::string CustomStructured::className() {
    return "atlas.grid.CustomStructured";
}


std::string CustomStructured::grid_type_str() {
    return "custom_structured";
}


CustomStructured::CustomStructured(const eckit::Parametrisation& params) :
    Structured(Domain::makeGlobal()) {
    setup(params);

    //if( ! params.get("short_name",shortName_) ) throw eckit::BadParameter("short_name missing in Params",Here());
    //if( ! params.has("hash") ) throw eckit::BadParameter("hash missing in Params",Here());
}


void CustomStructured::setup(const eckit::Parametrisation& params) {
    eckit::ValueList list;

    std::vector<long> pl;
    std::vector<double> latitudes;
    std::vector<double> lonmin;
    std::vector<double> lonmax;

    if( ! params.get("pl",pl) )
        throw eckit::BadParameter("pl missing in Params",Here());
    if( ! params.get("latitudes",latitudes) )
        throw eckit::BadParameter("latitudes missing in Params",Here());

    // Optionally specify N identifier
    params.get("N",Structured::N_);

    if (params.get("lon_min",lonmin) && params.get("lon_max",lonmax)) {
        Structured::setup(latitudes.size(),latitudes.data(),pl.data(),lonmin.data(),lonmax.data());
    }
    else {
        Structured::setup(latitudes.size(),latitudes.data(),pl.data());
    }
}


CustomStructured::CustomStructured(
    size_t nlat,
    const double lats[],
    const long nlons[],
    const Domain& dom ) :
    Structured(dom) {
    Structured::setup(nlat,lats,nlons);
}


CustomStructured::CustomStructured(
    size_t nlat,
    const double latitudes[],
    const long pl[],
    const double lonmin[],
    const double lonmax[],
    const Domain& dom ) :
    Structured(dom) {
    Structured::setup(nlat,latitudes,pl,lonmin,lonmax);
}


eckit::Properties CustomStructured::spec() const {
    eckit::Properties grid_spec;
    grid_spec.set("grid_type",  gridType());
    grid_spec.set("short_name", shortName());
    grid_spec.set("nlat",       nlat());
    grid_spec.set("latitudes",  eckit::makeVectorValue(latitudes()));
    grid_spec.set("pl",         eckit::makeVectorValue(pl()));
    grid_spec.set("lon_min",    eckit::makeVectorValue(lonmin_));
    grid_spec.set("lon_max",    eckit::makeVectorValue(lonmax_));
    if (N() != 0) {
        grid_spec.set("N", N());
    }
    return grid_spec;
}


extern "C" {


    Structured* atlas__grid__CustomStructured_int(size_t nlat, double lats[], int pl[]) {
        std::vector<long> pl_vector;
        pl_vector.assign(pl,pl+nlat);
        return new CustomStructured(nlat, lats, pl_vector.data(), Domain::makeGlobal());
    }


    Structured* atlas__grid__CustomStructured_long(size_t nlat, double lats[], long pl[]) {
        return new CustomStructured(nlat, lats, pl, Domain::makeGlobal());
    }


    Structured* atlas__grid__CustomStructured_lonmin_lonmax_int(size_t nlat, double lats[], int pl[], double lonmin[], double lonmax[]) {
        std::vector<long> pl_vector;
        pl_vector.assign(pl, pl+nlat);
        return new CustomStructured(nlat, lats, pl_vector.data(), lonmin, lonmax, Domain::makeGlobal());
    }


    Structured* atlas__grid__CustomStructured_lonmin_lonmax_long(size_t nlat, double lats[], long pl[], double lonmin[], double lonmax[]) {
        return new CustomStructured(nlat, lats, pl, lonmin, lonmax, Domain::makeGlobal());
    }


}


}  // namespace grid
}  // namespace atlas

