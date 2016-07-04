/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/gaussian/RegularGaussian.h"

#include "atlas/grid/gaussian/latitudes/Latitudes.h"


namespace atlas {
namespace grid {
namespace gaussian {


register_BuilderT1(Grid, RegularGaussian,RegularGaussian::grid_type_str());


std::string RegularGaussian::grid_type_str() {
    return "regular_gaussian";
}


std::string RegularGaussian::className() {
    return "atlas.grid.gaussian.RegularGaussian";
}


RegularGaussian::RegularGaussian(const eckit::Parametrisation& params) :
    Gaussian(Domain::makeGlobal()) {
    size_t N;
    if (!params.get("N",N)) {
        throw eckit::BadParameter("N missing in Params", Here());
    }
    setup(N,domain_);
}


RegularGaussian::RegularGaussian(const size_t& N, const Domain& dom) :
    Gaussian(dom) {
    setup(N,dom);
}


eckit::Properties RegularGaussian::spec() const {
    eckit::Properties grid_spec;
    grid_spec.set("grid_type", grid_type_str());
    grid_spec.set("short_name", shortName());
    grid_spec.set("N", N());
    grid_spec.set("nlat", nlat());
    return grid_spec;
}


void RegularGaussian::setup(const size_t& N, const Domain& dom) {

    // set internal Gaussian N and latitudes (assuming global domain)
    ASSERT(N>=2);
    Structured::N_  = N;
    Structured::lat_.assign(2*N, 0.);
    latitudes::gaussian_latitudes_npole_spole(N, Structured::lat_.data());

    // set (Ni,Nj) specific to domain-bound regular Gaussian
    const double lat_middle = (dom.north() + dom.south())/2.;
    const double lon_middle = (dom.east()  + dom.west() )/2.;
    const double inc_west_east = 90.0/static_cast<double>(N);

    Ni_ = N*4;
    if (!dom.isPeriodicEastWest()) {
        Ni_ = 0;
        for (size_t i=0; i<N*4; ++i) {
            const double lon = dom.west() + static_cast<double>(i*90.0)/static_cast<double>(N);
            if (dom.contains(lon,lat_middle))
                ++Ni_;
        }
    }
    ASSERT(0<Ni_ && Ni_<=N*4);

    Nj_ = N*2;
    if (!dom.includesPoleNorth() || !dom.includesPoleSouth()) {
        std::vector<double> lat_global;
        Structured::lat_.swap(lat_global);
        Structured::lat_.reserve(2*N);
        for (size_t i=0; i<N*2; ++i) {
            if (dom.contains(lon_middle,lat_global[i]))
                Structured::lat_.push_back(lat_global[i]);
        }
        Nj_ = Structured::lat_.size();
    }
    ASSERT(0<Nj_ && Nj_<=N*2);


    // set Structured:: members
    Structured::pl_.assign(Nj_, Ni_);
    Structured::nlonmin_ = Ni_;
    Structured::nlonmax_ = Ni_;
    Structured::npts_    = Ni_*Nj_;
    Structured::lon_inc_.assign(Nj_, inc_west_east);
    Structured::lonmin_ .assign(Nj_, dom.west());
    Structured::lonmax_ .assign(Nj_, dom.east() - (dom.isPeriodicEastWest()? inc_west_east : 0.));


    set_typeinfo();
}


void RegularGaussian::set_typeinfo() {
    std::stringstream s;
    s << "F"<< N();
    shortName_ = s.str();
    grid_type_ = grid_type_str();
}


extern "C" {


    Structured* atlas__grid__gaussian__RegularGaussian(size_t N) {
        return new RegularGaussian(N, Domain::makeGlobal());
    }


}


}  // namespace gaussian
}  // namespace grid
}  // namespace atlas

