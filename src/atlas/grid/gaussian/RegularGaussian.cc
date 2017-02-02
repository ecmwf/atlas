/*
 * (C) Copyright 1996-2017 ECMWF.
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


eckit::ConcreteBuilderT1<Grid, RegularGaussian> builder_RegularGaussian(RegularGaussian::grid_type_str());


std::string RegularGaussian::grid_type_str() {
    return "regular_gaussian";
}


std::string RegularGaussian::className() {
    return "atlas.grid.gaussian.RegularGaussian";
}


std::string RegularGaussian::gridType() const {
    return grid_type_str();
}


std::string RegularGaussian::shortName() const {
    if (shortName_.empty()) {
        std::stringstream s;
        s << "F" << N();
        if (!domain_.isGlobal()) {
            s << "-local";
        }
        shortName_ = s.str();
    }
    return shortName_;
}


RegularGaussian::RegularGaussian(const size_t& N, const Domain& domain) :
    Gaussian() {
    domain_ = domain;
    setup(N, domain);
}


RegularGaussian::RegularGaussian(const eckit::Parametrisation& params) :
    Gaussian() {
    size_t N;
    if (!params.get("N",N)) {
        throw eckit::BadParameter("N missing in Params", Here());
    }

    domain_ = Domain::makeGlobal();
    std::vector<double> p_domain(4);
    if (params.get("domain", p_domain)) {
      domain_ = Domain(p_domain[0], p_domain[1], p_domain[2], p_domain[3]);
    }

    setup(N, domain_);
}



namespace {
static eckit::Value domain_spec(const Domain& dom)
{
  std::vector<double> dspec(4);
  dspec[0] = dom.north();
  dspec[1] = dom.west();
  dspec[2] = dom.south();
  dspec[3] = dom.east();
  return eckit::makeVectorValue(dspec);
}
}

eckit::Properties RegularGaussian::spec() const {
    eckit::Properties grid_spec;
    grid_spec.set("grid_type", grid_type_str());
    grid_spec.set("short_name", shortName());
    grid_spec.set("N", N());
    grid_spec.set("nlat", nlat());
    grid_spec.set("domain", domain_spec(domain_) );
    return grid_spec;
}


void RegularGaussian::setup(const size_t& N, const Domain& dom) {

    // set internal Gaussian N and latitudes (assuming global domain)
    ASSERT(N>=2);
    Structured::N_  = N;
    std::vector<double> glb_lats(2*N);
    LatitudesNorthPoleToSouthPole(N,glb_lats.data());

    // set (Ni,Nj) specific to domain-bound regular Gaussian
    const double lon_middle = (dom.east()  + dom.west() )/2.;
    const double inc_west_east = 90.0/static_cast<double>(N);

    Nj_ = N*2;
    std::vector<double> dom_lats; dom_lats.reserve(2*N);
    if (!dom.includesPoleNorth() || !dom.includesPoleSouth()) {
        for (size_t i=0; i<N*2; ++i) {
            if (dom.contains(lon_middle,glb_lats[i]))
                dom_lats.push_back(glb_lats[i]);
        }
        Nj_ = dom_lats.size();
    }
    else {
      dom_lats = glb_lats;
    }
    ASSERT(0<Nj_ && Nj_<=N*2);

    double lonmin =  std::numeric_limits<double>::max();
    double lonmax = -std::numeric_limits<double>::max();
    const double tol = 1.e-6;
    Ni_ = 0;
    for( long jlon=-long(4*N)-1; jlon < long(4*N); ++jlon )
    {
      double lon = inc_west_east*jlon;
      if( lon+tol > dom.west() && lon-tol < dom.east() )
      {
          ++Ni_;
          lonmin = std::min(lonmin,lon);
          lonmax = std::max(lonmax,lon);
      }
    }
    if( dom.isPeriodicEastWest() )
    {
      Ni_=4*N;
      lonmax = lonmin +(Ni_-1)*inc_west_east;
    }

    std::vector<long>   dom_pl(Nj_,Ni_);
    std::vector<double> dom_lonmin(Nj_,lonmin);
    std::vector<double> dom_lonmax(Nj_,lonmax);
    Structured::setup(Nj_,dom_lats.data(),dom_pl.data(),dom_lonmin.data(),dom_lonmax.data());
}


extern "C" {


    Structured* atlas__grid__gaussian__RegularGaussian(size_t N) {
        return new RegularGaussian(N, Domain::makeGlobal());
    }


}


}  // namespace gaussian
}  // namespace grid
}  // namespace atlas

