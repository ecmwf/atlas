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

const std::vector<double>& ReducedGrid::lat() const
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

} // namespace atlas

