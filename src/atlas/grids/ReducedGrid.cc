/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <typeinfo> // std::bad_cast
#include <eckit/memory/Builder.h>
#include <eckit/memory/Factory.h>

#include "atlas/grids/ReducedGrid.h"
#include "atlas/GridSpec.h"
#include "atlas/GridSpecParams.h"
#include "atlas/util/Debug.h"

using eckit::Factory;
using eckit::Params;
using eckit::ValueParams;
using eckit::BadParameter;

namespace atlas {
namespace grids {


ReducedGrid* ReducedGrid::create( const Params& p )
{
  ReducedGrid* grid = dynamic_cast<ReducedGrid*>( Grid::create(p) );
  if( !grid )
    throw BadParameter("Grid is not a reduced grid",Here());
  return grid;
}

ReducedGrid* ReducedGrid::create(const std::string& uid)
{
  ReducedGrid* grid = dynamic_cast<ReducedGrid*>( Grid::create(uid) );
  if( !grid )
    throw BadParameter("Grid "+uid+" is not a reduced grid",Here());
  return grid;
}

ReducedGrid* ReducedGrid::create(const GridSpec& g)
{
  ReducedGrid* grid = dynamic_cast<ReducedGrid*>( Grid::create(g) );
  if( !grid )
    throw BadParameter("Grid is not a reduced grid",Here());
  return grid;
}

std::string ReducedGrid::className()
{
  return "atlas.ReducedGrid";
}

ReducedGrid::ReducedGrid() : N_(0)
{
}

ReducedGrid::ReducedGrid(const Params& params) : N_(0)
{
  setup(params);
  mask(params);

  if( ! params.has("grid_type") ) throw BadParameter("grid_type missing in Params",Here());
  if( ! params.has("uid") ) throw BadParameter("uid missing in Params",Here());
  if( ! params.has("hash") ) throw BadParameter("hash missing in Params",Here());

  grid_type_ = params.get("grid_type").as<std::string>();
  uid_ = params.get("uid").as<std::string>();
  hash_ = params.get("hash").as<std::string>();
}

void ReducedGrid::setup(const eckit::Params& params)
{
  eckit::ValueList list;

  std::vector<int> npts_per_lat;
  std::vector<double> latitudes;

  if( ! params.has("npts_per_lat") ) throw BadParameter("npts_per_lat missing in Params",Here());
  if( ! params.has("latitudes") ) throw BadParameter("latitudes missing in Params",Here());

  list = params.get("npts_per_lat");
  npts_per_lat.resize( list.size() );
  for(int j=0; j<npts_per_lat.size(); ++j)
    npts_per_lat[j] = list[j];

  list = params.get("latitudes");
  latitudes.resize( list.size() );
  for(int j=0; j<latitudes.size(); ++j)
    latitudes[j] = list[j];

  if( params.has("N") )
    N_ = params["N"];

  setup(latitudes.size(),latitudes.data(),npts_per_lat.data());
}

ReducedGrid::ReducedGrid( const std::vector<double>& _lats, const std::vector<size_t>& _nlons )
{
  int nlat = _nlons.size();
  std::vector<int> nlons(nlat);
  for( int j=0; j<nlat; ++j )
  {
    nlons[j] = static_cast<int>(nlons[j]);
  }
  setup(nlat,_lats.data(),nlons.data());
}

ReducedGrid::ReducedGrid( const int nlat, const double lats[], const int nlons[] )
{
  setup(nlat,lats,nlons);
}


void ReducedGrid::setup( const int nlat, const double lats[], const int nlons[], const double lonmin[], const double lonmax[] )
{
  ASSERT( nlat > 1 ); // can't have a grid with just one latitude

  nlons_.assign(nlons,nlons+nlat);

  lat_.assign(lats,lats+nlat);

  lonmin_.assign(lonmin,lonmin+nlat);
  lonmax_.assign(lonmax,lonmax+nlat);

  npts_ = 0;
  nlonmax_ = 0;
  double lon_min(1000), lon_max(-1000);

  for( int jlat=0; jlat<nlat; ++jlat )
  {
    //ASSERT( nlon(jlat) > 1 ); // can't have grid with just one longitude
    nlonmax_ = std::max(nlon(jlat),nlonmax_);

    lon_min = std::min(lon_min,lonmin_[jlat]);
    lon_max = std::max(lon_max,lonmax_[jlat]);

    npts_ += nlons_[jlat];
  }

  // Default is global
  bounding_box_ = BoundBox(lat_[0]/*north*/, lat_[nlat-1]/*south*/, lon_max/*east*/, lon_min/*west*/ );
}


void ReducedGrid::setup( const int nlat, const double lats[], const int nlons[] )
{
  std::vector<double> lonmin(nlat,0.);
  std::vector<double> lonmax(nlat);
  for( int jlat=0; jlat<nlat; ++jlat )
    lonmax[jlat] = 360.-360./static_cast<double>(nlons[jlat]);
  setup(nlat,lats,nlons,lonmin.data(),lonmax.data());
}

void ReducedGrid::setup_colat_hemisphere(const int N, const double colat[], const int lon[], const AngleUnit unit)
{
  std::vector<int> nlons(2*N);
  std::copy( lon, lon+N, nlons.begin() );
  std::reverse_copy( lon, lon+N, nlons.begin()+N );
  std::vector<double> lats(2*N);
  colat_to_lat_hemisphere(N,colat,lats.data(),unit);
  std::reverse_copy( lats.data(), lats.data()+N, lats.begin()+N );
  double convert = (unit == RAD ? 180.*M_1_PI : 1.);
  for( int j=0; j<N; ++j )
    lats[j] *= convert;
  for( int j=N; j<2*N; ++j )
    lats[j] *= -convert;
  setup(2*N,lats.data(),nlons.data());
}

void ReducedGrid::setup_lat_hemisphere(const int N, const double lat[], const int lon[], const AngleUnit unit)
{
  std::vector<int> nlons(2*N);
  std::copy( lon, lon+N, nlons.begin() );
  std::reverse_copy( lon, lon+N, nlons.begin()+N );
  std::vector<double> lats(2*N);
  std::copy( lat, lat+N, lats.begin() );
  std::reverse_copy( lat, lat+N, lats.begin()+N );
  double convert = (unit == RAD ? 180.*M_1_PI : 1.);
  for( int j=0; j<N; ++j )
    lats[j] *= convert;
  for( int j=N; j<2*N; ++j )
    lats[j] *= -convert;
  setup(2*N,lats.data(),nlons.data());
}

int ReducedGrid::N() const
{
  if( N_==0 )
  {
    throw eckit::Exception("N cannot be returned because grid of type "+grid_type()+
                           " is not based on a global grid.", Here() );
  }
  return N_;
}

Grid::BoundBox ReducedGrid::bounding_box() const
{
  return bounding_box_;
}

size_t ReducedGrid::npts() const
{
  return npts_;
}

void ReducedGrid::lonlat( double crds[] ) const
{
  int c(0);
  for( int jlat=0; jlat<nlat(); ++jlat )
  {
    double ylat = lat(jlat);
    for( int jlon=0; jlon<nlon(jlat); ++jlon )
    {
      crds[c++] = lon(jlat,jlon);
      crds[c++] = ylat;
    }
  }
}

void ReducedGrid::lonlat( std::vector<Point>& pts ) const
{
  pts.resize(npts());
  int c(0);
  for( int jlat=0; jlat<nlat(); ++jlat )
  {
    double y = lat(jlat);
    for( int jlon=0; jlon<nlon(jlat); ++jlon )
    {
      pts[c++].assign(lon(jlat,jlon),y);
    }
  }
}

std::string ReducedGrid::grid_type() const
{
  return grid_type_;
}

GridSpec ReducedGrid::spec() const
{
  GridSpec grid_spec(grid_type());

  grid_spec.uid( uid() );
  grid_spec.set("hash", hash());
  grid_spec.set("nlat",nlat());
  grid_spec.set_latitudes(latitudes());
  grid_spec.set_npts_per_lat(npts_per_lat());
  grid_spec.set_bounding_box(bounding_box());

  if( N_ != 0 )
    grid_spec.set("N", N_ );

  return grid_spec;
}

bool ReducedGrid::same(const Grid& _g) const
{
  try
  {
    const ReducedGrid& g = *dynamic_cast<const ReducedGrid*>(&_g);

    if( npts() != g.npts() )
      return false;

    if( nlat() != g.nlat() )
      return false;

    for( int jlat=0; jlat<nlat(); ++jlat )
    {
      if( nlon(jlat) != g.nlon(jlat) )
        return false;
    }

    return true;
  }
  catch (std::bad_cast& e)
  {
    return false;
  }
  return false;
}

int ReducedGrid::nlat() const
{
  return lat_.size();
}

int ReducedGrid::nlon(int jlat) const
{
  return nlons_[jlat];
}

int ReducedGrid::nlonmax() const
{
  return nlonmax_;
}

const std::vector<int>&  ReducedGrid::npts_per_lat() const
{
  return nlons_;
}

double ReducedGrid::lon(const int jlat, const int jlon) const
{
  return lonmin_[jlat] + (double)jlon * (lonmax_[jlat]-lonmin_[jlat]) / ( (double)nlon(jlat) - 1. );
}

double ReducedGrid::lat(const int jlat) const
{
  return lat_[jlat];
}

void ReducedGrid::lonlat( const int jlon, const int jlat, double crd[] ) const
{
  crd[0] = lon(jlat,jlon);
  crd[1] = lat(jlat);
}

const std::vector<double>& ReducedGrid::latitudes() const
{
  return lat_;
}

std::string ReducedGrid::uid() const
{
  ASSERT( !uid_.empty() );
  return uid_;
}

std::string ReducedGrid::hash() const
{
  ASSERT( !hash_.empty() );
  return hash_;
}

void ReducedGrid::mask(const Domain& dom)
{
  // If the mask is larger or equal to the domain, no need to mask!
  if ( dom.max().lat() >= bounding_box_.max().lat() &&
       dom.max().lon() >= bounding_box_.max().lon() &&
       dom.min().lat() <= bounding_box_.min().lat() &&
       dom.min().lon() <= bounding_box_.min().lon() )
    return;

  /// 1) Figure out mask
  const double bblatmin = dom.min().lat()-degrees_eps();
  const double bblatmax = dom.max().lat()+degrees_eps();
  const double bblonmin = dom.min().lon()-degrees_eps();
  const double bblonmax = dom.max().lon()+degrees_eps();

  std::vector<int> masklat(nlat(),-1);
  std::vector<int> masklon_first(nlat(),-1);
  std::vector<int> masklon_last (nlat(),-1);
  int new_nlat(0);
  double new_bblatmin( 90.);
  double new_bblatmax(-90.);
  double new_bblonmin(360.);
  double new_bblonmax(  0.);
  const int __nlat = nlat();
  for( int jlat=0; jlat<__nlat; ++jlat )
  {
    masklat[jlat] = 0;
    const double __lat = lat(jlat);
    if( __lat >= bblatmin && __lat <= bblatmax )
    {
      const int __nlon = nlon(jlat);
      masklon_first[jlat] = __nlon;
      masklon_last [jlat] = -1;

//      if( __nlon == 0 )
//      {
//        if( !masklat[jlat] ) ++new_nlat;
//        masklat[jlat] = 1;
//        masklon_first[jlat] = 0;
//        masklon_last[jlat]  = -1;
//        new_bblatmin = std::min(new_bblatmin,__lat);
//        new_bblatmax = std::max(new_bblatmax,__lat);
//      }

      for( int jlon=0; jlon<__nlon; ++jlon )
      {
        const double __lon = lon(jlat,jlon);
        if( __lon >= bblonmin && __lon <= bblonmax )
        {
          if( !masklat[jlat] ) ++new_nlat;
          masklat[jlat] = 1;
          masklon_first[jlat] = std::min(masklon_first[jlat],jlon);
          masklon_last [jlat] = std::max(masklon_last [jlat],jlon);

          new_bblatmin = std::min(new_bblatmin,__lat);
          new_bblatmax = std::max(new_bblatmax,__lat);
          new_bblonmin = std::min(new_bblonmin,__lon);
          new_bblonmax = std::max(new_bblonmax,__lon);
        }
      }
    }
  }

  /// 2) modify grid
  std::vector<double> new_lat;     new_lat.reserve(new_nlat);
  std::vector<double> new_lonmin;  new_lonmin.reserve(new_nlat);
  std::vector<double> new_lonmax;  new_lonmax.reserve(new_nlat);
  std::vector<int>    new_nlons;   new_nlons.reserve(new_nlat);
  int new_npts(0);
  for( int jlat=0; jlat<__nlat; ++jlat )
  {
    if( masklat[jlat] )
    {
      new_lat.   push_back( lat(jlat) );
      new_lonmin.push_back( lon(jlat,masklon_first[jlat]) );
      new_lonmax.push_back( lon(jlat,masklon_last [jlat]) );
      new_nlons. push_back( masklon_last[jlat]-masklon_first[jlat]+1 );
      new_npts += new_nlons.back();
    }
  }

  npts_   = new_npts;
  lat_    = new_lat;
  nlons_  = new_nlons;
  lonmin_ = new_lonmin;
  lonmax_ = new_lonmax;
  bounding_box_ = BoundBox(new_bblatmax,new_bblatmin,new_bblonmax,new_bblonmin);
}

void ReducedGrid::mask(const eckit::Params& p)
{
  if( p.has("domain_s") )
  {
    Domain dom(p["domain_n"],p["domain_s"],p["domain_e"],p["domain_w"]);
    mask( dom );
  }
  else
  {
    // If there is a "bbox_s" value in params,
    // then make_bounding_box(params).global() will be false.
    mask( make_bounding_box(p) );
  }
}


int  atlas__ReducedGrid__nlat(ReducedGrid* This)
{
  return This->nlat();
}

void atlas__ReducedGrid__nlon(ReducedGrid* This, const int* &nlons, int &size)
{
  nlons = This->npts_per_lat().data();
  size  = This->npts_per_lat().size();
}

int atlas__ReducedGrid__npts(ReducedGrid* This)
{
  return This->npts();
}

double atlas__ReducedGrid__lon(ReducedGrid* This,int jlat,int jlon)
{
  return This->lon(jlat,jlon);
}

double atlas__ReducedGrid__lat(ReducedGrid* This,int jlat)
{
  return This->lat(jlat);
}

void atlas__ReducedGrid__latitudes(ReducedGrid* This, const double* &lat, int &size)
{
  lat  = This->latitudes().data();
  size = This->latitudes().size();
}

} // namespace grids
} // namespace atlas
