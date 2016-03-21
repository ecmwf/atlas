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

Structured::Structured() : Global(), N_(0)
{
}

//------------------------------------------------------------------------------

void Structured::setup(const eckit::Parametrisation& params)
{
  eckit::ValueList list;

  std::vector<long> npts_per_lat;
  std::vector<double> latitudes;

  if( ! params.get("npts_per_lat",npts_per_lat) ) throw BadParameter("npts_per_lat missing in Params",Here());
  if( ! params.get("latitudes",latitudes) ) throw BadParameter("latitudes missing in Params",Here());

  params.get("N",N_);

  setup(latitudes.size(),latitudes.data(),npts_per_lat.data());
}

//------------------------------------------------------------------------------

void Structured::setup( const size_t nlat, const double lats[], const long nlons[], const double lonmin[], const double lonmax[] )
{
  ASSERT(nlat > 1);  // can't have a grid with just one latitude

  nlons_.assign(nlons,nlons+nlat);

  lat_.assign(lats,lats+nlat);

  lonmin_.assign(lonmin,lonmin+nlat);
  lonmax_.assign(lonmax,lonmax+nlat);

  lon_inc_.resize(nlat);
  
  npts_ = 0;
  nlonmax_ = 0;
  nlonmin_ = std::numeric_limits<size_t>::max();
  double lon_min(1000), lon_max(-1000);

  for(size_t jlat = 0; jlat < nlat; ++jlat)
  {
    //ASSERT( nlon(jlat) > 1 ); // can't have grid with just one longitude
    nlonmax_ = std::max(nlon(jlat),nlonmax_);
    nlonmin_ = std::min(nlon(jlat),nlonmin_);

    lon_min = std::min(lon_min,lonmin_[jlat]);
    lon_max = std::max(lon_max,lonmax_[jlat]);
    lon_inc_[jlat] = (lonmax_[jlat] - lonmin_[jlat])/(nlons_[jlat]-1);

    npts_ += nlons_[jlat];
  }

  bounding_box_ = BoundBox(lat_[0]/*north*/, lat_[nlat-1]/*south*/, lon_max/*east*/, lon_min/*west*/ );
}

//------------------------------------------------------------------------------

void Structured::setup( const size_t nlat, const double lats[], const long nlons[] )
{
  std::vector<double> lonmin(nlat,0.);
  std::vector<double> lonmax(nlat);
  for(size_t jlat = 0; jlat < nlat; ++jlat)
  {
    if( nlons[jlat] )
      lonmax[jlat] = 360.-360./static_cast<double>(nlons[jlat]);
    else
      lonmax[jlat] = 0.;
  }
  setup(nlat,lats,nlons,lonmin.data(),lonmax.data());
}

//------------------------------------------------------------------------------

void Structured::setup_lat_hemisphere(const size_t N, const double lat[], const long lon[])
{
  std::vector<long> nlons(2*N);
  std::copy( lon, lon+N, nlons.begin() );
  std::reverse_copy( lon, lon+N, nlons.begin()+N );
  std::vector<double> lats(2*N);
  std::copy( lat, lat+N, lats.begin() );
  std::reverse_copy( lat, lat+N, lats.begin()+N );
  for(size_t j = N; j < 2*N; ++j)
    lats[j] *= -1.;
  setup(2*N,lats.data(),nlons.data());
}

//------------------------------------------------------------------------------

size_t Structured::N() const
{
  return N_;
}

//------------------------------------------------------------------------------

BoundBox Structured::boundingBox() const
{
  return bounding_box_;
}

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
  bounding_box_.hash(md5);
}

//------------------------------------------------------------------------------

extern "C" {

size_t atlas__grid__global__Structured__N(Structured* This)
{
  return This->N();
}

size_t atlas__grid__global__Structured__nlat(Structured* This)
{
  return This->nlat();
}

size_t atlas__grid__global__Structured__nlon(Structured* This, size_t jlat)
{
  return This->nlon(jlat);
}

void atlas__grid__global__Structured__pl(Structured* This, const long* &nlons, size_t &size)
{
  nlons = This->pl().data();
  size  = This->pl().size();
}

size_t atlas__grid__global__Structured__nlonmax(Structured* This)
{
  return This->nlonmax();
}

size_t atlas__grid__global__Structured__nlonmin(Structured* This)
{
  return This->nlonmin();
}

size_t atlas__grid__global__Structured__npts(Structured* This)
{
  return This->npts();
}

double atlas__grid__global__Structured__lat(Structured* This,size_t jlat)
{
  return This->lat(jlat);
}

double atlas__grid__global__Structured__lon(Structured* This,size_t jlat,size_t jlon)
{
  return This->lon(jlat, jlon);
}

void atlas__grid__global__Structured__lonlat(Structured* This, size_t jlat, size_t jlon, double crd[])
{
  This->lonlat(jlat, jlon, crd);
}

void atlas__grid__global__Structured__latitudes(Structured* This, const double* &lat, size_t &size)
{
  lat  = This->latitudes().data();
  size = This->latitudes().size();
}

int atlas__grid__global__Structured__reduced  (Structured* This)
{
  return This->reduced();
}


Structured* atlas__grid__global__Structured(char* identifier)
{
  return Structured::create( std::string(identifier) );
}

void atlas__grid__global__Structured__delete(Structured* This)
{
  delete This;
}

}

//------------------------------------------------------------------------------

} // namespace global
} // namespace grid
} // namespace atlas
