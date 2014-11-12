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

#define DEBUG_OUTPUT 0

namespace atlas {
namespace meshgen {

static const double rad_to_deg = 180.*M_1_PI;

void colat_to_lat_hemisphere(const int N, const double colat[], double lats[], const AngleUnit unit)
{
  std::copy( colat, colat+N, lats );
  double pole = (unit == DEG ? 90. : M_PI_2);
  for (int i=0; i<N; ++i)
    lats[i]=pole-lats[i];
}

void predict_gaussian_colatitudes_hemisphere(const int N, double colat[])
{
  double z;
  for( int i=0; i<N; ++i )
  {
    z = (4.*(i+1.)-1.)*M_PI/(4.*2.*N+2.);
    colat[i] = ( z+1./(tan(z)*(8.*(2.*N)*(2.*N))) ) * rad_to_deg;
  }
}

void predict_gaussian_latitudes_hemisphere(const int N, double lats[])
{
  std::vector<double> colat(N);
  predict_gaussian_colatitudes_hemisphere(N,colat.data());
  colat_to_lat_hemisphere(N,colat.data(),lats,DEG);
}


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

void RGG::setup_colat_hemisphere(const int N, const int lon[], const double colat[], const AngleUnit unit)
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

void RGG::setup_lat_hemisphere(const int N, const int lon[], const double lat[], const AngleUnit unit)
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

GG::GG(int nlon, int N) : RGG()
{
  std::vector<double> lat(N);
  std::vector<int>    nlons(N,nlon);
  predict_gaussian_latitudes_hemisphere(N,lat.data());
  setup_lat_hemisphere(N,nlons.data(),lat.data(),DEG);
}


RegularGrid::RegularGrid(int nlon, int nlat) : RGG()
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

RGG* new_reduced_gaussian_grid( const std::string& identifier )
{
  RGG* rgg = 0;

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

RGG* new_reduced_gaussian_grid( const std::vector<long>& nlon )
{
  return new RGG(nlon.size(),nlon.data());
}

// ------------------------------------------------------------------

RGG* new_regular_latlon_grid( int nlon, int nlat )
{
  if( nlon%2 != 0 ) throw eckit::BadParameter("nlon must be even number",Here());
  if( nlat%2 != 0 ) throw eckit::BadParameter("nlat must be even number",Here());
  return new RegularGrid(nlon,nlat);
}

// ------------------------------------------------------------------

RGG* new_regular_gaussian_grid( int nlon, int nlat )
{
  if( nlon%2 != 0 ) throw eckit::BadParameter("nlon must be even number",Here());
  if( nlat%2 != 0 ) throw eckit::BadParameter("nlat must be even number",Here());
  return new GG(nlon,nlat/2.);
}


RGG* atlas__new_reduced_gaussian_grid(char* identifier)
{
	return new_reduced_gaussian_grid( std::string(identifier) );
}

RGG* atlas__new_regular_gaussian_grid ( int nlon, int nlat )
{
  return new_regular_gaussian_grid( nlon, nlat );
}

RGG* atlas__new_regular_latlon_grid(int nlon, int nlat)
{
  return new_regular_latlon_grid( nlon, nlat );
}

RGG* atlas__new_custom_reduced_gaussian_grid(int nlon[], int nlat)
{
  std::vector<long> nlon_vector;
  nlon_vector.assign(nlon,nlon+nlat);
  return new_reduced_gaussian_grid(nlon_vector);
}

int  atlas__RGG__nlat(RGG* This)
{
  return This->nlat();
}

void atlas__RGG__nlon(RGG* This, const int* &nlons, int &size)
{
  nlons = This->nlons().data();
  size  = This->nlons().size();
}

int atlas__RGG__ngptot(RGG* This)
{
  return This->nPoints();
}

double atlas__RGG__lon(RGG* This,int jlon,int jlat)
{
  return This->lon(jlat,jlon);
}

double atlas__RGG__lat(RGG* This,int jlat)
{
  return This->lat(jlat);
}

void atlas__RGG__lats(RGG* This, const double* &lat, int &size)
{
  lat  = This->lat().data();
  size = This->lat().size();
}


} // namespace meshgen
} // namespace atlas

