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

RGG::RGG(const int nlat, const int lon[])
{
  /*
  First prediction of colatitudes
  */
  std::vector<double> colat(nlat);
  double z;
  for( int i=0; i<nlat; ++i )
  {
    z = (4.*(i+1.)-1.)*M_PI/(4.*2.*nlat+2.);
    colat[i] = z+1./(tan(z)*(8.*(2.*nlat)*(2.*nlat)));
  }
  /*
  Fill in final structures
  */
  lat_.resize(2*nlat);
  lon_.resize(2*nlat);
  std::copy( lon, lon+nlat, lon_.begin() );
  std::reverse_copy( lon, lon+nlat, lon_.begin()+nlat );
  std::copy( colat.begin(), colat.begin()+nlat, lat_.begin() );
  std::reverse_copy( colat.begin(), colat.begin()+nlat, lat_.begin()+nlat );
  for (int i=0; i<nlat; ++i)
    lat_[i]=M_PI/2.-lat_[i];
  for (int i=nlat; i<2*nlat; ++i)
    lat_[i]=-M_PI/2.+lat_[i];
}
RGG::RGG(const size_t nlat, const long lon[])
{
  /*
  First prediction of colatitudes
  */
  std::vector<double> colat(nlat);
  double z;
  for( int i=0; i<nlat; ++i )
  {
    z = (4.*(i+1.)-1.)*M_PI/(4.*2.*nlat+2.);
    colat[i] = z+1./(tan(z)*(8.*(2.*nlat)*(2.*nlat)));
  }
  /*
  Fill in final structures
  */
  lat_.resize(2*nlat);
  lon_.resize(2*nlat);
  std::copy( lon, lon+nlat, lon_.begin() );
  std::reverse_copy( lon, lon+nlat, lon_.begin()+nlat );
  std::copy( colat.begin(), colat.begin()+nlat, lat_.begin() );
  std::reverse_copy( colat.begin(), colat.begin()+nlat, lat_.begin()+nlat );
  for (int i=0; i<nlat; ++i)
    lat_[i]=M_PI/2.-lat_[i];
  for (int i=nlat; i<2*nlat; ++i)
    lat_[i]=-M_PI/2.+lat_[i];
}


int RGG::ngptot() const
{
  return std::accumulate(lon_.data(),lon_.data()+lon_.size(),0);
}

GG::GG(int nlon, int nlat) : RGG()
{
  /*
  First prediction of colatitudes
  */
  std::vector<double> colat(nlat);
  double z;
  for( int i=0; i<nlat; ++i )
  {
    z = (4.*(i+1.)-1.)*M_PI/(4.*2.*nlat+2.);
    colat[i] = z+1./(tan(z)*(8.*(2.*nlat)*(2.*nlat)));
  }
  /*
  Fill in final structures
  */
  lon_.assign(2*nlat,nlon);
  lat_.resize(2*nlat);
  std::copy( colat.begin(), colat.begin()+nlat, lat_.begin() );
  std::reverse_copy( colat.begin(), colat.begin()+nlat, lat_.begin()+nlat );
  for (int i=0; i<nlat; ++i)
    lat_[i]=M_PI/2.-lat_[i];
  for (int i=nlat; i<2*nlat; ++i)
    lat_[i]=-M_PI/2.+lat_[i];
}


RegularGrid::RegularGrid(int nlon, int nlat) : RGG()
{
  double dy = M_PI/static_cast<double>(nlat);

  lon_.resize(nlat);
  lon_.assign(nlat,nlon);

  lat_.resize(nlat);
  for( int i=0; i<nlat; ++i )
  {
    lat_[i] = M_PI_2 - (i+0.5)*dy;
  }
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

void atlas__RGG__nlon(RGG* This, const int* &nlon, int &size)
{
  nlon = This->nlon().data();
  size = This->nlon().size();
}

int atlas__RGG__ngptot(RGG* This)
{
  return This->ngptot();
}

double atlas__RGG__lon(RGG* This,int jlon,int jlat)
{
  return This->lon(jlon,jlat);
}

double atlas__RGG__lat(RGG* This,int jlat)
{
  return This->lat(jlat);
}

void atlas__RGG__lats(RGG* This, const double* &lats, int &size)
{
  lats = This->lat().data();
  size = This->lat().size();
}


} // namespace meshgen
} // namespace atlas

