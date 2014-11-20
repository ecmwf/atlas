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

#include "atlas/ReducedGrid.h"
#include "atlas/GridSpec.h"
#include "atlas/util/Debug.h"

namespace atlas {

namespace {
static double units =1;// M_PI/180.; // radians
}

std::string ReducedGrid::className()
{
  return "atlas.ReducedGrid";
}

ReducedGrid::ReducedGrid()
{
}

ReducedGrid::ReducedGrid(const eckit::Params& params)
{
  eckit::ValueList list;

  int N;
  std::vector<int> npts_per_lat;
  std::vector<double> latitudes;

  if( ! params.has("npts_per_lat") ) throw eckit::BadParameter("npts_per_lat missing in Params",Here());
  if( ! params.has("latitudes") ) throw eckit::BadParameter("latitudes missing in Params",Here());
  if( ! params.has("N") ) throw eckit::BadParameter("N missing in Params",Here());
  if( ! params.has("grid_type") ) throw eckit::BadParameter("grid_type missing in Params",Here());
  if( ! params.has("uid") ) throw eckit::BadParameter("uid missing in Params",Here());
  if( ! params.has("hash") ) throw eckit::BadParameter("hash missing in Params",Here());

  list = params.get("npts_per_lat");
  npts_per_lat.resize( list.size() );
  for(int j=0; j<npts_per_lat.size(); ++j)
    npts_per_lat[j] = list[j];

  list = params.get("latitudes");
  latitudes.resize( list.size() );
  for(int j=0; j<latitudes.size(); ++j)
    latitudes[j] = list[j];

  setup(latitudes.size(),npts_per_lat.data(),latitudes.data());

  grid_type_ = params.get("grid_type").as<std::string>();
  uid_ = params.get("uid").as<std::string>();
  hash_ = params.get("hash").as<std::string>();
}

ReducedGrid::ReducedGrid( const std::vector<size_t>& _nlons, const std::vector<double>& _lats )
{
  int nlat = _nlons.size();
  std::vector<int> nlons(nlat);
  for( int j=0; j<nlat; ++j )
  {
    nlons[j] = static_cast<int>(nlons[j]);
  }
  setup(nlat,nlons.data(),_lats.data());
}

ReducedGrid::ReducedGrid( const int nlons[], const double lats[], const int nlat )
{
  setup(nlat,nlons,lats);
}

void ReducedGrid::setup(const int nlat, const int nlons[], const double lats[])
{
  N_ = nlat/2;
  nlons_.assign(nlons,nlons+nlat);

  lat_.assign(lats,lats+nlat);

  npts_ = 0;
  nlonmax_ = 0;
  double lon_min(1000), lon_max(-1000);
  for( int jlat=0; jlat<nlat; ++jlat )
  {
    nlonmax_ = std::max(nlon(jlat),nlonmax_);
    lon_min = std::min(lon_min,lon(jlat,0));
    lon_max = std::max(lon_max,lon(jlat,nlon(jlat)-1));
    npts_ += nlons_[jlat];
  }

  bounding_box_ = BoundBox(lat_[0]/*north*/, lat_[nlat-1]/*south*/, lon_max/*east*/, lon_min/*west*/ );

  if( grid_type_.empty() )
    grid_type_ = "ReducedGrid";
}

void ReducedGrid::setup_colat_hemisphere(const int N, const int lon[], const double colat[], const AngleUnit unit)
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
  setup(2*N,nlons.data(),lats.data());
}

void ReducedGrid::setup_lat_hemisphere(const int N, const int lon[], const double lat[], const AngleUnit unit)
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
  setup(2*N,nlons.data(),lats.data());
}


Grid::BoundBox ReducedGrid::bounding_box() const
{
  return bounding_box_;
}

size_t ReducedGrid::npts() const
{
  return npts_;
}

void ReducedGrid::coordinates( std::vector<double>& crds) const
{
  ///@todo this should be first longitude, then latitude
  crds.resize(2*npts());
  int c(0);
  for( int jlat=0; jlat<nlat(); ++jlat )
  {
    double y = lat(jlat);
    for( int jlon=0; jlon<nlon(jlat); ++jlat )
    {
      crds[c++] = y;
      crds[c++] = lon(jlat,jlon);
    }
  }
}

void ReducedGrid::coordinates( std::vector<Point>& pts ) const
{
  ///@todo this should be first longitude, then latitude
  pts.resize(npts());
  int c(0);
  for( int jlat=0; jlat<nlat(); ++jlat )
  {
    double y = lat(jlat);
    for( int jlon=0; jlon<nlon(jlat); ++jlat )
    {
      pts[c++].assign(y,lon(jlat,jlon));
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
  grid_spec.set("hash", hash_);
  grid_spec.set("N", N_ );
  grid_spec.set_bounding_box(bounding_box_);
  grid_spec.set_latitudes(latitudes());
  grid_spec.set_npts_per_lat(npts_per_lat());

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
  return static_cast<double>(jlon) * 360./static_cast<double>(nlon(jlat)) * units;
}

double ReducedGrid::lat(const int jlat) const
{
  return lat_[jlat] * units;
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

double atlas__ReducedGrid__lon(ReducedGrid* This,int jlon,int jlat)
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

} // namespace atlas

