/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/Structured.h"

#include <algorithm>
#include <limits>
#include "atlas/runtime/ErrorHandling.h"


namespace atlas {
namespace grid {


Structured* Structured::create(const eckit::Parametrisation& params) {
    Structured* grid = dynamic_cast<Structured*>(Grid::create(params));
    if (!grid)
        throw eckit::BadParameter("Grid is not a reduced grid", Here());
    return grid;

}


Structured* Structured::create(const std::string& uid) {
    Structured* grid = dynamic_cast<Structured*>( Grid::create(uid) );
    if (!grid)
        throw eckit::BadParameter("Grid "+uid+" is not a reduced grid",Here());
    return grid;
}


std::string Structured::className() {
    return "atlas.grid.Structured";
}


std::string Structured::grid_type_str() {
    return "structured";
}


Structured::Structured() :
    Grid(),
    N_(0) {
}


Structured::~Structured() {
}


void Structured::setup(
    const size_t nlat,
    const double lats[],
    const long pl[],
    const double lonmin[],
    const double lonmax[] ) {
    ASSERT(nlat>1);  // can't have a grid with just one latitude

    pl_    .assign(pl, pl+nlat);
    lat_   .assign(lats, lats+nlat);
    lonmin_.assign(lonmin, lonmin+nlat);
    lonmax_.assign(lonmax, lonmax+nlat);
    npts_ = static_cast<size_t>(std::accumulate(pl_.begin(), pl_.end(), 0));

    lon_inc_.resize(nlat);
    nlonmin_ = nlonmax_ = static_cast<size_t>(pl_[0]);

    for (size_t jlat = 0; jlat < nlat; ++jlat) {
        lon_inc_[jlat] = (lonmax_[jlat]-lonmin_[jlat])/double(pl_[jlat]-1);
        nlonmin_ = std::min(static_cast<size_t>(pl_[jlat]),nlonmin_);
        nlonmax_ = std::max(static_cast<size_t>(pl_[jlat]),nlonmax_);
    }
}


void Structured::setup_lon_limits(const size_t nlat, const long pl[], const Domain& dom, double lonmin[], double lonmax[]) {
    ASSERT(nlat);

    std::fill_n(lonmin, nlat, dom.west());
    std::fill_n(lonmax, nlat, dom.east());

    const double ew = dom.east() - dom.west();
    const bool isPeriodicEastWest = dom.isPeriodicEastWest();

    for (size_t jlat = 0; jlat < nlat; ++jlat) {
        const double ndiv = static_cast<double>(pl[jlat] + (isPeriodicEastWest? 0:-1));
        lonmax[jlat] -= pl[jlat]? ew/ndiv : 0.;;
    }
}


void Structured::setup_lat_hemisphere(const size_t N, const double lat[], const long lon[]) {

    const size_t nlat = 2*N;
    ASSERT(nlat);

    // construct pl (assuming global domain)
    std::vector<long> pl(nlat);
    std::copy( lon, lon+N, pl.begin() );
    std::reverse_copy( lon, lon+N, pl.begin()+static_cast<long>(N) );

    // construct latitudes (assuming global domain)
    std::vector<double> lats(nlat);
    std::copy( lat, lat+N, lats.begin() );
    std::reverse_copy( lat, lat+N, lats.begin()+static_cast<long>(N) );
    for(size_t j = N; j < nlat; ++j) {
        lats[j] *= -1.;
    }

    // assign longitude limits
    std::vector<double> lonmin(nlat);
    std::vector<double> lonmax(nlat);
    setup_lon_limits(nlat, pl.data(), domain(), lonmin.data(), lonmax.data());

    // common setup
    setup(nlat, lats.data(), pl.data(), lonmin.data(), lonmax.data());
}


void Structured::lonlat( std::vector<Point>& pts ) const {
    pts.resize(npts());

    for(size_t jlat=0, c=0; jlat<nlat(); ++jlat) {
        const double y = lat(jlat);
        for(size_t jlon=0; jlon<nlon(jlat); ++jlon) {
            pts[c++].assign(lon(jlat,jlon),y);
        }
    }
}


size_t Structured::copyLonLatMemory(double* pts, size_t size) const {
    size_t sizePts = 2*npts();
    ASSERT(size >= sizePts);

    for(size_t jlat=0, c=0; jlat<nlat(); ++jlat ) {
        const double y = lat(jlat);
        for( size_t jlon=0; jlon<nlon(jlat); ++jlon ) {
            pts[c++] = lon(jlat,jlon);
            pts[c++] = y;
        }
    }
    return sizePts;
}


void Structured::print(std::ostream& os) const {
    os << "Structured(Name:" << shortName() << ")";
}


void Structured::hash(eckit::MD5& md5) const {
    // Through inheritance the grid_type_str() might differ while still being same grid
    //md5.add(grid_type_str());

    md5.add(latitudes().data(), sizeof(double)*latitudes().size());
    md5.add(pl().data(), sizeof(long)*nlat());
}


// --------------------------------------------------------------------


extern "C" {


    size_t atlas__grid__Structured__N(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->N();
        );
        return 0;
    }


    size_t atlas__grid__Structured__nlat(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nlat();
        );
        return 0;
    }


    size_t atlas__grid__Structured__nlon(Structured* This, size_t jlat) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nlon(jlat);
        );
        return 0;
    }


    void atlas__grid__Structured__pl(Structured* This, const long* &nlons, size_t &size) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            nlons = This->pl().data();
            size  = This->pl().size();
        );
    }


    size_t atlas__grid__Structured__nlonmax(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nlonmax();
        );
        return 0;
    }


    size_t atlas__grid__Structured__nlonmin(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->nlonmin();
        );
        return 0;
    }


    size_t atlas__grid__Structured__npts(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->npts();
        );
        return 0;
    }


    double atlas__grid__Structured__lat(Structured* This,size_t jlat) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->lat(jlat);
        );
        return 0.;
    }


    double atlas__grid__Structured__lon(Structured* This,size_t jlat,size_t jlon) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->lon(jlat, jlon);
        );
        return 0.;
    }


    void atlas__grid__Structured__lonlat(Structured* This, size_t jlat, size_t jlon, double crd[]) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            This->lonlat(jlat, jlon, crd);
        );
    }


    void atlas__grid__Structured__latitudes(Structured* This, const double* &lat, size_t &size) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            lat  = This->latitudes().data();
            size = This->latitudes().size();
        );
    }


    int atlas__grid__Structured__reduced(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
            return This->reduced();
        );
        return 1;
    }


    Structured* atlas__grid__Structured(char* identifier) {
        ATLAS_ERROR_HANDLING(
            ASSERT( identifier );
            return Structured::create( std::string(identifier) );
        );
        return 0;
    }


    Structured* atlas__grid__Structured__config(eckit::Parametrisation* conf) {
        ATLAS_ERROR_HANDLING(
            ASSERT( conf );
            return Structured::create(*conf);
        );
        return 0;
    }


    void atlas__grid__Structured__delete(Structured* This) {
        ATLAS_ERROR_HANDLING(
            ASSERT( This );
        );
        delete This;
    }


}


}  // namespace grid
}  // namespace atlas

