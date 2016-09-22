/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/gaussian/Gaussian.h"

#include <typeinfo>
#include "eckit/memory/Builder.h"
#include "atlas/grid/gaussian/latitudes/Latitudes.h"

namespace atlas {
namespace grid {
namespace gaussian {


void Gaussian::LatitudesNorthPoleToEquator(const size_t N, double latitudes[]) {
    gaussian::latitudes::gaussian_latitudes_npole_equator(N,latitudes);
}


void Gaussian::LatitudesNorthPoleToSouthPole(const size_t N, double latitudes[]) {
    gaussian::latitudes::gaussian_latitudes_npole_spole(N,latitudes);
}


void Gaussian::QuadratureNorthPoleToEquator(const size_t N, double weights[]) {
    std::vector<double> lats(N);
    latitudes::gaussian_quadrature_npole_equator(N,lats.data(),weights);
}


void Gaussian::QuadratureNorthPoleToSouthPole(const size_t N, double weights[]) {
    std::vector<double> lats(2*N);
    latitudes::gaussian_quadrature_npole_equator(N,lats.data(),weights);
}


std::string Gaussian::grid_type_str() {
    return "gaussian";
}


std::string Gaussian::className() {
    return "atlas.grid.gaussian.Gaussian";
}


Gaussian::Gaussian() :
    Structured() {
}


void Gaussian::setup_N_hemisphere(const size_t N, const long pl[], const domain::Domain& dom) {
    Structured::N_ = N;

    if( dom.isGlobal() && dom.west() == 0. ) {
        // hemisphere
        std::vector<double> lats (N);
        LatitudesNorthPoleToEquator(N,lats.data());
        Structured::setup_lat_hemisphere(N,lats.data(),pl);
    }
    else
    {
        std::vector<double> glb_lats (2*N);
        std::vector<long>   glb_pl   (2*N);
        for( size_t jlat=0; jlat<N; ++jlat )
        {
            glb_pl[jlat] = pl[jlat];
            glb_pl[2*N-1-jlat] = pl[jlat];
        }
        LatitudesNorthPoleToSouthPole(N,glb_lats.data());
        std::vector<double> dom_lats;   dom_lats.  reserve(2*N);
        std::vector<long>   dom_pl;     dom_pl.    reserve(2*N);
        std::vector<double> dom_lonmin; dom_lonmin.reserve(2*N);
        std::vector<double> dom_lonmax; dom_lonmax.reserve(2*N);
        const double tol = 1.e-6;
        size_t nlat = 0;
        const bool periodic_east_west = dom.isPeriodicEastWest();
        for( size_t jlat=0; jlat<2*N; ++jlat )
        {
            if( glb_lats[jlat]-tol < dom.north() && glb_lats[jlat]+tol > dom.south() )
            {
                ++nlat;
                const double lat = glb_lats[jlat];
                double lonmin =  std::numeric_limits<double>::max();
                double lonmax = -std::numeric_limits<double>::max();
                size_t nlon = 0;
                const double inc_west_east = 360./double(glb_pl[jlat]);
                for( long jlon=-glb_pl[jlat]; jlon<glb_pl[jlat]; ++jlon )
                {
                    const double lon = inc_west_east*jlon;
                    if( lon+tol > dom.west() && lon-tol < dom.east() )
                    {
                        ++nlon;
                        lonmin = std::min(lonmin,lon);
                        lonmax = std::max(lonmax,lon);
                    }
                }
                if( periodic_east_west )
                {
                  nlon = glb_pl[jlat];
                  lonmax = lonmin + (nlon-1)*inc_west_east;
                }
                dom_lats  .push_back(lat);
                dom_pl    .push_back(nlon);
                dom_lonmin.push_back(lonmin);
                dom_lonmax.push_back(lonmax);
            }
        }
        Structured::setup(nlat,dom_lats.data(),dom_pl.data(),dom_lonmin.data(),dom_lonmax.data());
    }
}


eckit::Properties Gaussian::spec() const {
    eckit::Properties grid_spec;

    grid_spec.set("grid_type", gridType() );
    grid_spec.set("short_name",shortName());

    grid_spec.set("nlat",nlat());
    grid_spec.set("N", N() );

    grid_spec.set("pl",eckit::makeVectorValue(pl()));

    return grid_spec;
}


}  // namespace gaussian
}  // namespace grid
}  // namespace atlas
