/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <typeinfo>
#include "eckit/memory/Builder.h"
#include "atlas/internals/atlas_config.h"
#include "atlas/grid/global/lonlat/ShiftedLon.h"

using eckit::BadParameter;
using eckit::Params;

namespace atlas {
namespace grid {
namespace global {
namespace lonlat {

//------------------------------------------------------------------------------

register_BuilderT1(Grid,ShiftedLon,ShiftedLon::grid_type_str());

//------------------------------------------------------------------------------

std::string ShiftedLon::grid_type_str()
{
  return "shifted_lon";
}

//------------------------------------------------------------------------------

std::string ShiftedLon::className()
{
  return "atlas.grid.global.lonlat.ShiftedLon";
}

//------------------------------------------------------------------------------

void ShiftedLon::set_typeinfo()
{
  std::stringstream s;
  if( N() ) {
    s << "Slon" << N();
  } else {
    s << "Slon" << nlon() << "x" << nlat();
  }
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

//------------------------------------------------------------------------------

ShiftedLon::ShiftedLon(const eckit::Parametrisation& p)
  : LonLat(Shift::LON)
{
  setup(p);
  set_typeinfo();
}

//------------------------------------------------------------------------------

ShiftedLon::ShiftedLon( const size_t N )
  : LonLat(Shift::LON)
{
  setup(N);
  set_typeinfo();
}

//------------------------------------------------------------------------------

ShiftedLon::ShiftedLon( const int nlon, const int nlat )
  : LonLat(Shift::LON)
{
  setup( (size_t)nlon, (size_t)nlat );
  set_typeinfo();
}

//------------------------------------------------------------------------------

ShiftedLon::ShiftedLon( const size_t nlon, const size_t nlat )
  : LonLat(Shift::LON)
{
  setup( (size_t)nlon, (size_t) nlat);
  set_typeinfo();
}

//------------------------------------------------------------------------------

ShiftedLon::ShiftedLon( const double &londeg, const double &latdeg )
  : LonLat(Shift::LON)
{
  setup(londeg,latdeg);
  set_typeinfo();
}

//------------------------------------------------------------------------------

void ShiftedLon::setup(const eckit::Parametrisation& p)
{
  size_t nlon, nlat, N(0);
  p.get("N",N);

  if( N > 0 )
  {
    setup(N);
  }
  else
  {
    if( !p.has("nlon") && !p.has("lon_inc") ) throw BadParameter("nlon or lon_inc missing in Params",Here());
    if( !p.has("nlat") && !p.has("lat_inc") ) throw BadParameter("nlat or lat_inc missing in Params",Here());

    double lon_inc, lat_inc;
    if (p.get("nlon",nlon) && p.get("nlat",nlat))
    {
      setup(nlon,nlat);
    }
    else if (p.get("lon_inc",lon_inc) && p.get("lat_inc",lat_inc))
    {
      setup(lon_inc,lat_inc);
    }
    else
    {
      throw BadParameter("Bad combination of parameters");
    }
  }}

//------------------------------------------------------------------------------

void ShiftedLon::setup( const size_t N )
{
  double delta = 90./static_cast<double>(N);
  std::vector<double> lats(2*N+1);
  std::vector<long>   nlons(2*N+1,4*N);
  std::vector<double> lonmin(2*N+1,0.+0.5*delta);
  std::vector<double> lonmax(2*N+1,360.-0.5*delta);

  double latmax = 90.;

  for( size_t jlat=0; jlat<2*N+1; ++jlat )
  {
    lats[jlat] = latmax - static_cast<double>(jlat)*delta;
  }

  Structured::N_ = N;
  Structured::setup(2*N+1,lats.data(),nlons.data(),lonmin.data());
}

//------------------------------------------------------------------------------

void ShiftedLon::setup( const size_t nlon, const size_t nlat )
{
  double latmin = -90.;
  double latmax = +90.;
  double londeg = 360./static_cast<double>(nlon);
  double latdeg = (latmax-latmin)/static_cast<double>(nlat-1);

  std::vector<double> lats(nlat);
  std::vector<long>   nlons(nlat,nlon);
  std::vector<double> lonmin(nlat,0.5*londeg);
  std::vector<double> lonmax(nlat,360.-0.5*londeg);

  for( size_t jlat=0; jlat<nlat; ++jlat )
  {
    lats[jlat] = latmax - static_cast<double>(jlat)*latdeg;
  }

  if( (nlat-1)%2 == 0 && nlon==2*(nlat-1) )
  {
    Structured::N_ = (nlat-1)/2;
  }
  Structured::setup(nlat,lats.data(),nlons.data(),lonmin.data());
}

//------------------------------------------------------------------------------

void ShiftedLon::setup( const double londeg, const double latdeg )
{
  double Llon = 360.-londeg;
  double Llat = 180.;
  double nlon_real = Llon/londeg + 1.;
  double nlat_real = Llat/latdeg + 1.;
  size_t nlon = static_cast<size_t>(nlon_real);
  size_t nlat = static_cast<size_t>(nlat_real);
  if( nlon_real - nlon > 0. )
  {
    std::stringstream msg;
    msg << Llon << " is not divisible by londeg " << londeg << " --> nlon = " << nlon_real;
    throw BadParameter(msg.str(),Here());
  }
  if( nlat_real - nlat > 0. )
  {
    std::stringstream msg;
    msg << Llat << " is not divisible by latdeg " << latdeg << " --> nlat = " << nlat_real;
    throw BadParameter(msg.str(),Here());
  }
  setup(nlon,nlat);
}

//------------------------------------------------------------------------------

eckit::Properties ShiftedLon::spec() const
{
  eckit::Properties grid_spec;

  grid_spec.set("grid_type",gridType() );
  grid_spec.set("short_name",shortName());

  grid_spec.set("N", N() );
  grid_spec.set("nlon", nlon() );
  grid_spec.set("nlat", nlat() );

  BoundBox bbox = boundingBox();
  grid_spec.set("bbox_s", bbox.min().lat());
  grid_spec.set("bbox_w", bbox.min().lon());
  grid_spec.set("bbox_n", bbox.max().lat());
  grid_spec.set("bbox_e", bbox.max().lon());

  return grid_spec;
}

//------------------------------------------------------------------------------

extern "C"
{

Structured* atlas__grid__global__lonlat__ShiftedLon(size_t nlon, size_t nlat)
{
  return new ShiftedLon(nlon,nlat);
}

}

//-----------------------------------------------------------------------------

} // namespace global
} // namespace lonlat
} // namespace grid
} // namespace atlas
