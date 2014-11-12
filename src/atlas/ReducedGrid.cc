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
#include "atlas/ReducedGrid.h"
#include "atlas/GridSpec.h"

namespace atlas {

namespace {
static double units = M_PI/180.; // radians
}

std::string ReducedGrid::className()
{
  return "atlas.ReducedGrid";
}

ReducedGrid::ReducedGrid()
{
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
  bbox_ = BoundBox(lat_[0]/*north*/, lat_[nlat-1]/*south*/, lon_max/*east*/, lon_min/*west*/ );
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


Grid::BoundBox ReducedGrid::boundingBox() const
{
  return bbox_;
}

size_t ReducedGrid::nPoints() const
{
  return npts_;
}

void ReducedGrid::coordinates( std::vector<double>& crds) const
{
  ///@todo this should be first longitude, then latitude
  crds.resize(2*nPoints());
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
  pts.resize(nPoints());
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

std::string ReducedGrid::gridType() const
{
  return ReducedGrid::className();
}

GridSpec ReducedGrid::spec() const
{
  return GridSpec(ReducedGrid::className());
}

bool ReducedGrid::same(const Grid& _g) const
{
  try
  {
    const ReducedGrid& g = *dynamic_cast<const ReducedGrid*>(&_g);

    if( nPoints() != g.nPoints() )
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

const std::vector<int>&  ReducedGrid::nlons() const
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
  throw eckit::NotImplemented("atlas::ReducedGrid is a base class and has no unique id",Here());
  return ReducedGrid::className();
}

std::string ReducedGrid::hash() const
{
  throw eckit::NotImplemented("atlas::ReducedGrid is a base class and has no hash",Here());
  return ReducedGrid::className();
}

int  atlas__ReducedGrid__nlat(ReducedGrid* This)
{
  return This->nlat();
}

void atlas__ReducedGrid__nlon(ReducedGrid* This, const int* &nlons, int &size)
{
  nlons = This->nlons().data();
  size  = This->nlons().size();
}

int atlas__ReducedGrid__npts(ReducedGrid* This)
{
  return This->nPoints();
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

