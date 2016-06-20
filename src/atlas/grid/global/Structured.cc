/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <limits>
#include <typeinfo>
#include <string>
#include "eckit/memory/Builder.h"
#include "eckit/memory/Factory.h"
#include "atlas/grid/global/Structured.h"
#include "atlas/internals/Debug.h"
#include "atlas/runtime/ErrorHandling.h"

using eckit::Factory;
using eckit::MD5;
using eckit::BadParameter;

namespace atlas {
namespace grid {
namespace global {

//------------------------------------------------------------------------------

Structured* Structured::create(const eckit::Parametrisation& p) {

  Structured* grid = dynamic_cast<Structured*>(Grid::create(p));
  if (!grid) throw BadParameter("Grid is not a reduced grid", Here());
  return grid;

}

//------------------------------------------------------------------------------

Structured* Structured::create(const std::string& uid)
{
  Structured* grid = dynamic_cast<Structured*>( Grid::create(uid) );
  if( !grid )
    throw BadParameter("Grid "+uid+" is not a reduced grid",Here());
  return grid;
}

//------------------------------------------------------------------------------

std::string Structured::className() { return "atlas.grid.global.Structured"; }

//------------------------------------------------------------------------------

Structured::Structured(const Domain& dom) :
    Global(dom),
    N_(0)
{
}

//------------------------------------------------------------------------------

void Structured::setup(const size_t nlat, const double lats[], const long pl[], const double lonmin[] )
{
  ASSERT(nlat > 1);  // can't have a grid with just one latitude

  pl_.assign(pl,pl+nlat);

  lat_.assign(lats,lats+nlat);

  lonmin_.assign(lonmin,lonmin+nlat);
  lon_inc_.resize(nlat);

  npts_    = 0;
  nlonmax_ = 0;
  nlonmin_ = std::numeric_limits<size_t>::max();
  double lon_min(1000), lon_max(-1000);

  const double ew = domain_.east()  - domain_.west() ;
  const bool isPeriodicEastWest = domain_.isPeriodicEastWest();

  for(size_t jlat = 0; jlat < nlat; ++jlat)
  {
    const double ndiv = static_cast<double>(pl_[jlat] + (isPeriodicEastWest? 0:-1));
    lon_inc_[jlat] = pl_[jlat]? ew/ndiv : 0.;


    lon_min = std::min(lon_min,lonmin_[jlat]);
    lon_max = std::max(lon_max,lonmin_[jlat]+(pl_[jlat]-1l)*lon_inc_[jlat]);
    nlonmax_ = std::max( (size_t) pl_[jlat],nlonmax_);
    nlonmin_ = std::min( (size_t) pl_[jlat],nlonmin_);

    npts_ += pl_[jlat];
  }

}

//------------------------------------------------------------------------------

void Structured::setup( const size_t nlat, const double lats[], const long pl[] )
{
  std::vector<double> lonmin(nlat,0.);
  setup(nlat,lats,pl,lonmin.data());
}

//------------------------------------------------------------------------------

void Structured::setup_lat_hemisphere(const size_t N, const double lat[], const long lon[])
{
  std::vector<long> pl(2*N);
  std::copy( lon, lon+N, pl.begin() );
  std::reverse_copy( lon, lon+N, pl.begin()+N );
  std::vector<double> lats(2*N);
  std::copy( lat, lat+N, lats.begin() );
  std::reverse_copy( lat, lat+N, lats.begin()+N );
  for(size_t j = N; j < 2*N; ++j)
    lats[j] *= -1.;
  setup(2*N,lats.data(),pl.data());
}

//------------------------------------------------------------------------------

size_t Structured::N() const
{
  return N_;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

size_t Structured::npts() const { return npts_; }

//------------------------------------------------------------------------------

void Structured::lonlat( std::vector<Point>& pts ) const
{
  pts.resize(npts());
  int c(0);
  for(size_t jlat = 0; jlat < nlat(); ++jlat)
  {
    double y = lat(jlat);
    for(size_t jlon = 0; jlon < nlon(jlat); ++jlon)
    {
      pts[c++].assign(lon(jlat,jlon),y);
    }
  }
}

//------------------------------------------------------------------------------

std::string Structured::gridType() const
{
  return grid_type_;
}

//------------------------------------------------------------------------------

std::string Structured::getOptimalMeshGenerator() const
{
    return "Structured";
}

//------------------------------------------------------------------------------

size_t Structured::copyLonLatMemory(double* pts, size_t size) const
{
    size_t sizePts = 2*npts();

    ASSERT(size >= sizePts);

    for(size_t c = 0, jlat=0; jlat<nlat(); ++jlat )
    {
      double y = lat(jlat);
      for( size_t jlon=0; jlon<nlon(jlat); ++jlon )
      {
        pts[c++] = lon(jlat,jlon);
        pts[c++] = y;
      }
    }
    return sizePts;
}

//------------------------------------------------------------------------------

void Structured::print(std::ostream& os) const
{
    os << "Structured(Name:" << shortName() << ")";
}

//------------------------------------------------------------------------------

const std::vector<double>& Structured::latitudes() const
{
  return lat_;
}

//------------------------------------------------------------------------------

std::string Structured::shortName() const {
  ASSERT(!shortName_.empty());
  return shortName_;
}

//------------------------------------------------------------------------------

void Structured::hash(eckit::MD5& md5) const {
  // Through inheritance the grid_type_str() might differ while still being same grid
      //md5.add(grid_type_str());

  md5.add(latitudes().data(), sizeof(double)*latitudes().size());
  md5.add(pl().data(), sizeof(long)*nlat());
  domain_.hash(md5);
}

//------------------------------------------------------------------------------

extern "C" {

size_t atlas__grid__global__Structured__N(Structured* This)
{
  ATLAS_ERROR_HANDLING(
     ASSERT( This );
     return This->N();
  );
  return 0;
}

size_t atlas__grid__global__Structured__nlat(Structured* This)
{
  ATLAS_ERROR_HANDLING(
     ASSERT( This );
     return This->nlat();
  );
  return 0;
}

size_t atlas__grid__global__Structured__nlon(Structured* This, size_t jlat)
{
  ATLAS_ERROR_HANDLING(
     ASSERT( This );
     return This->nlon(jlat);
  );
  return 0;
}

void atlas__grid__global__Structured__pl(Structured* This, const long* &nlons, size_t &size)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    nlons = This->pl().data();
    size  = This->pl().size();
  );
}

size_t atlas__grid__global__Structured__nlonmax(Structured* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->nlonmax();
  );
  return 0;
}

size_t atlas__grid__global__Structured__nlonmin(Structured* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->nlonmin();
  );
  return 0;
}

size_t atlas__grid__global__Structured__npts(Structured* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->npts();
  );
  return 0;
}

double atlas__grid__global__Structured__lat(Structured* This,size_t jlat)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->lat(jlat);
  );
  return 0.;
}

double atlas__grid__global__Structured__lon(Structured* This,size_t jlat,size_t jlon)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->lon(jlat, jlon);
  );
  return 0.;
}

void atlas__grid__global__Structured__lonlat(Structured* This, size_t jlat, size_t jlon, double crd[])
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    This->lonlat(jlat, jlon, crd);
  );
}

void atlas__grid__global__Structured__latitudes(Structured* This, const double* &lat, size_t &size)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    lat  = This->latitudes().data();
    size = This->latitudes().size();
  );
}

int atlas__grid__global__Structured__reduced  (Structured* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->reduced();
  );
  return 1;
}


Structured* atlas__grid__global__Structured(char* identifier)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( identifier );
    return Structured::create( std::string(identifier) );
  );
  return 0;
}

Structured* atlas__grid__global__Structured__config(eckit::Parametrisation* conf)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( conf );
    return Structured::create(*conf);
  );
  return 0;
}

void atlas__grid__global__Structured__delete(Structured* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
  );
  delete This;
}

}

//------------------------------------------------------------------------------

} // namespace global
} // namespace grid
} // namespace atlas
