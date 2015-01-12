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
#include "eckit/memory/Factory.h"
#include "eckit/memory/Builder.h"

#include "atlas/grids/GaussianLatitudes.h"
#include "atlas/meshgen/RGG.h"
#include "atlas/Util.h"

namespace atlas {
namespace meshgen {

RGG::RGG(const int N, const int lon[]) : grids::ReducedGrid()
{
  std::vector<double> lat(N);
  grids::gaussian_latitudes_npole_equator(N,lat.data());
  setup_lat_hemisphere(N,lat.data(),lon,DEG);
}

RGG::RGG(const size_t N, const long lon[])
{
  std::vector<double> lat(N);
  std::vector<int>    nlons(N);
  std::copy( lon, lon+N, nlons.begin() );
  grids::gaussian_latitudes_npole_equator(N,lat.data());
  setup_lat_hemisphere(N,lat.data(),nlons.data(),DEG);
}


GG::GG(int nlon, int N)
{
  std::vector<double> lat(N);
  std::vector<int>    nlons(N,nlon);
  grids::gaussian_latitudes_npole_equator(N,lat.data());
  setup_lat_hemisphere(N,lat.data(),nlons.data(),DEG);
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
  setup(nlat,lats.data(),nlons.data());
}

// ------------------------------------------------------------------

grids::ReducedGrid* new_reduced_gaussian_grid( const std::string& identifier )
{
  grids::ReducedGrid* rgg = 0;

  if     ( identifier == "T63"   ) rgg = new T63();
  else if( identifier == "T95"   ) rgg = new T95();
  else if( identifier == "T159"  ) rgg = new T159();
  else if( identifier == "T255"  ) rgg = new T255();
  else if( identifier == "T511"  ) rgg = new T511();
  else if( identifier == "T1279" ) rgg = new T1279();
  else if( identifier == "T2047" ) rgg = new T2047();
  else if( identifier == "T3999" ) rgg = new T3999();
  else if( identifier == "T7999" ) rgg = new T7999();
  else
  {
    if( !eckit::Factory<Grid>::instance().exists(identifier) )
    {
      std::stringstream msg;
      msg << "Cannot find grid "<<identifier<<" in " << eckit::Factory<grids::ReducedGrid>::instance();
      throw eckit::BadParameter(msg.str(),Here());
    }
    rgg = grids::ReducedGrid::create(identifier);
  }
  return rgg;
}

// ------------------------------------------------------------------

grids::ReducedGrid* new_reduced_gaussian_grid( const std::vector<long>& nlon )
{
  return new RGG(nlon.size(),nlon.data());
}

// ------------------------------------------------------------------

grids::ReducedGrid* new_regular_latlon_grid( int nlon, int nlat )
{
  if( nlon%2 != 0 ) throw eckit::BadParameter("nlon must be even number",Here());
  if( nlat%2 != 0 ) throw eckit::BadParameter("nlat must be even number",Here());
  return new RegularGrid(nlon,nlat);
}

// ------------------------------------------------------------------

grids::ReducedGrid* new_regular_gaussian_grid( int nlon, int nlat )
{
  if( nlon%2 != 0 ) throw eckit::BadParameter("nlon must be even number",Here());
  if( nlat%2 != 0 ) throw eckit::BadParameter("nlat must be even number",Here());
  return new GG(nlon,nlat/2.);
}


grids::ReducedGrid* atlas__new_reduced_gaussian_grid(char* identifier)
{
	return new_reduced_gaussian_grid( std::string(identifier) );
}

grids::ReducedGrid* atlas__new_regular_gaussian_grid ( int nlon, int nlat )
{
  return new_regular_gaussian_grid( nlon, nlat );
}

grids::ReducedGrid* atlas__new_regular_latlon_grid(int nlon, int nlat)
{
  return new_regular_latlon_grid( nlon, nlat );
}

grids::ReducedGrid* atlas__new_custom_reduced_gaussian_grid(int nlon[], int nlat)
{
  std::vector<long> nlon_vector;
  nlon_vector.assign(nlon,nlon+nlat);
  return new_reduced_gaussian_grid(nlon_vector);
}

} // namespace meshgen
} // namespace atlas

