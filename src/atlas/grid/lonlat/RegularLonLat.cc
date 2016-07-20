/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/lonlat/RegularLonLat.h"


namespace atlas {
namespace grid {
namespace lonlat {


register_BuilderT1(Grid,RegularLonLat,RegularLonLat::grid_type_str());


std::string RegularLonLat::grid_type_str() {
    return "regular_lonlat";
}


std::string RegularLonLat::className() {
    return "atlas.grid.lonlat.RegularLonLat";
}


void RegularLonLat::set_typeinfo() {
    std::stringstream s;
    if( N() ) {
        s << "L" << N();
    } else {
        s << "L" << nlon() << "x" << nlat();
    }
    shortName_ = s.str();
    grid_type_ = grid_type_str();
}


RegularLonLat::RegularLonLat(const eckit::Parametrisation& p) :
    LonLat(Shift::NONE,Domain::makeGlobal()) {
    setup(p);
}


RegularLonLat::RegularLonLat(const size_t N, const Domain& dom) :
    LonLat(Shift::NONE, dom) {
    LonLat::setup(N, dom);
}


RegularLonLat::RegularLonLat(const size_t nlon, const size_t nlat, const Domain& dom) :
    LonLat(Shift::NONE, dom) {
    LonLat::setup(nlon, nlat, dom);
}


void RegularLonLat::setup(const eckit::Parametrisation& p) {
    size_t nlon, nlat, N(0);

    std::vector<double> p_domain(4);

    if( p.get("domain", p_domain) )
    {
      domain_ = Domain(p_domain[0],p_domain[1],p_domain[2],p_domain[3]);
    }
    else
    {
      domain_ = Domain::makeGlobal();
    }

    if (p.get("N",N)) {
        LonLat::setup(N, domain_);
    }
    else if (p.get("nlon",nlon) && p.get("nlat",nlat)) {
        LonLat::setup(nlon, nlat, domain_);
    }
    else {
        throw eckit::BadParameter("Params (nlon,nlat) or N missing",Here());
    }
}


eckit::Properties RegularLonLat::spec() const {
    eckit::Properties grid_spec;
    grid_spec.set("grid_type",  gridType());
    grid_spec.set("short_name", shortName());
    grid_spec.set("N",    N());
    grid_spec.set("nlon", nlon());
    grid_spec.set("nlat", nlat());
    grid_spec.set("domain", domain_spec(domain_) );
    return grid_spec;
}


extern "C" {


    Structured* atlas__grid__lonlat__RegularLonLat(size_t nlon, size_t nlat) {
        return new RegularLonLat(nlon,nlat);
    }


}


}  // namespace lonlat
}  // namespace grid
}  // namespace atlas

