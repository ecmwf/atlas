/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <numeric>
#include <algorithm>
#include <eckit/exception/Exceptions.h>
#include "atlas/meshgen/RGG.h"
#include "atlas/Util.h"

#define DEBUG_OUTPUT 0

namespace atlas {
namespace meshgen {

RGG::RGG(const int N, const int lon[]) : ReducedGrid()
{
  std::vector<double> lat(N);
  predict_gaussian_latitudes_hemisphere(N,lat.data());
  setup_lat_hemisphere(N,lon,lat.data(),DEG);
}

RGG::RGG(const size_t N, const long lon[])
{
  std::vector<double> lat(N);
  std::vector<int>    nlons(N);
  std::copy( lon, lon+N, nlons.begin() );
  predict_gaussian_latitudes_hemisphere(N,lat.data());
  setup_lat_hemisphere(N,nlons.data(),lat.data(),DEG);

}


GG::GG(int nlon, int N)
{
  std::vector<double> lat(N);
  std::vector<int>    nlons(N,nlon);
  predict_gaussian_latitudes_hemisphere(N,lat.data());
  setup_lat_hemisphere(N,nlons.data(),lat.data(),DEG);
}


RegularGrid::RegularGrid(int nlon, int nlat)
{
  std::vector<double> lats(nlat);
  std::vector<int>    nlons(nlat,nlon);

  double dy = 180./static_cast<double>(nlat);
  for( int i=0; i<nlat; ++i )
  {
    lats[i] = 90. - (i+0.5)*dy;
  }
  setup(nlat,nlons.data(),lats.data());
}

// ------------------------------------------------------------------

ReducedGrid* new_reduced_gaussian_grid( const std::string& identifier )
{
  ReducedGrid* rgg = 0;

  if     ( identifier == "T63"   ) rgg = new T63();
  else if( identifier == "T95"   ) rgg = new T95();
  else if( identifier == "T159"  ) rgg = new T159();
  else if( identifier == "T255"  ) rgg = new T255();
  else if( identifier == "T511"  ) rgg = new T511();
  else if( identifier == "T1279" ) rgg = new T1279();
  else if( identifier == "T2047" ) rgg = new T2047();
  else if( identifier == "T3999" ) rgg = new T3999();
  else if( identifier == "T7999" ) rgg = new T7999();
  else throw eckit::BadParameter("Cannot create grid"+identifier,Here());
  return rgg;
}

// ------------------------------------------------------------------

ReducedGrid* new_reduced_gaussian_grid( const std::vector<long>& nlon )
{
  return new RGG(nlon.size(),nlon.data());
}

// ------------------------------------------------------------------

ReducedGrid* new_regular_latlon_grid( int nlon, int nlat )
{
  if( nlon%2 != 0 ) throw eckit::BadParameter("nlon must be even number",Here());
  if( nlat%2 != 0 ) throw eckit::BadParameter("nlat must be even number",Here());
  return new RegularGrid(nlon,nlat);
}

// ------------------------------------------------------------------

ReducedGrid* new_regular_gaussian_grid( int nlon, int nlat )
{
  if( nlon%2 != 0 ) throw eckit::BadParameter("nlon must be even number",Here());
  if( nlat%2 != 0 ) throw eckit::BadParameter("nlat must be even number",Here());
  return new GG(nlon,nlat/2.);
}


ReducedGrid* atlas__new_reduced_gaussian_grid(char* identifier)
{
	return new_reduced_gaussian_grid( std::string(identifier) );
}

ReducedGrid* atlas__new_regular_gaussian_grid ( int nlon, int nlat )
{
  return new_regular_gaussian_grid( nlon, nlat );
}

ReducedGrid* atlas__new_regular_latlon_grid(int nlon, int nlat)
{
  return new_regular_latlon_grid( nlon, nlat );
}

ReducedGrid* atlas__new_custom_reduced_gaussian_grid(int nlon[], int nlat)
{
  std::vector<long> nlon_vector;
  nlon_vector.assign(nlon,nlon+nlat);
  return new_reduced_gaussian_grid(nlon_vector);
}

int  atlas__RGG__nlat(ReducedGrid* This)
{
  return This->nlat();
}

void atlas__RGG__nlon(ReducedGrid* This, const int* &nlons, int &size)
{
  nlons = This->nlons().data();
  size  = This->nlons().size();
}

int atlas__RGG__ngptot(ReducedGrid* This)
{
  return This->nPoints();
}

double atlas__RGG__lon(ReducedGrid* This,int jlon,int jlat)
{
  return This->lon(jlat,jlon);
}

double atlas__RGG__lat(ReducedGrid* This,int jlat)
{
  return This->lat(jlat);
}

void atlas__RGG__lats(ReducedGrid* This, const double* &lat, int &size)
{
  lat  = This->lat().data();
  size = This->lat().size();
}


} // namespace meshgen
} // namespace atlas

