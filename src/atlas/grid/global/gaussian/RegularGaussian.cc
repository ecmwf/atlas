/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/global/gaussian/RegularGaussian.h"

#include <typeinfo>
#include "eckit/memory/Builder.h"
#include "eckit/types/FloatCompare.h"
#include "atlas/grid/global/gaussian/latitudes/Latitudes.h"


namespace atlas {
namespace grid {
namespace global {
namespace gaussian {


register_BuilderT1(Grid, RegularGaussian,RegularGaussian::grid_type_str());


void RegularGaussian::set_typeinfo() {
    std::stringstream s;
    s << "F"<< N();
    shortName_ = s.str();
    grid_type_ = grid_type_str();
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

    grid_spec.set("grid_type",grid_type_str());
    grid_spec.set("short_name",shortName());

    grid_spec.set("N", N() );
    grid_spec.set("nlat",nlat());

    return grid_spec;
}


void RegularGaussian::setup(const size_t& N, const Domain& dom) {

    // set internal Gaussian N
    ASSERT(N>=2);
    N_ = N;

    // set (Ni,Nj) specific to domain-bound regular Gaussian
    Ni_ = N*4;
    if (!dom.isPeriodicEastWest()) {
        const double lat = dom.north();

        Ni_ = 0;
        for (size_t i=0; i<N*4; ++i) {
            const double lon = dom.west() + static_cast<double>(i*90.0)/static_cast<double>(N);
            if (dom.contains(lon,lat))
                ++ Ni_;
        }
    }
    ASSERT(0<Ni_ && Ni_<N*4);

    Nj_ = N*2;
    if (!dom.includesPoleNorth() || !dom.includesPoleSouth()) {
        const double lon = dom.west();

        std::vector<double> lats(2*N,0.);
        latitudes::gaussian_latitudes_npole_spole(N,lats.data());

        Nj_ = 0;
        for (size_t i=0; i<N*2; ++i) {
            if (dom.contains(lon,lats[i]))
                ++ Nj_;
        }
    }
    ASSERT(0<Nj_ && Nj_<N*2);

    // set internal Structured/Gaussian
    std::vector<long> pl(Nj_,Ni_);
    Gaussian::setup_N_hemisphere(N,pl.data());

    set_typeinfo();
}


extern "C" {
Structured* atlas__grid__global__gaussian__RegularGaussian(size_t N)
{
  return new RegularGaussian(N);
}
}


} // namespace gaussian
} // namespace global
} // namespace grid
} // namespace atlas
