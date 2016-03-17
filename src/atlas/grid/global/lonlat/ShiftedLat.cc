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
#include "atlas/grid/global/lonlat/ShiftedLat.h"

using eckit::BadParameter;
using eckit::Params;

namespace atlas {
namespace grid {
namespace global {
namespace lonlat {

//------------------------------------------------------------------------------

register_BuilderT1(Grid,ShiftedLat,ShiftedLat::grid_type_str());

//------------------------------------------------------------------------------

std::string ShiftedLat::grid_type_str()
{
  return "shifted_lat";
}

//------------------------------------------------------------------------------

std::string ShiftedLat::className()
{
  return "atlas.grid.global.lonlat.ShiftedLat";
}

//------------------------------------------------------------------------------

void ShiftedLat::set_typeinfo()
{
  std::stringstream s;
  if( N() ) {
    s << "Slat" << N();
  } else {
    s << "Slat" << nlon() << "x" << nlat();
  }

  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

//------------------------------------------------------------------------------

ShiftedLat::ShiftedLat(const eckit::Parametrisation& p) : LonLat()
{
  setup(p);
  set_typeinfo();
}

//------------------------------------------------------------------------------

ShiftedLat::ShiftedLat( const long N ) : LonLat()
{
  setup(N);
  set_typeinfo();
}

//------------------------------------------------------------------------------

ShiftedLat::ShiftedLat( const long nlon, const long nlat ) : LonLat()
{
  setup(nlon,nlat);
  set_typeinfo();
}

//------------------------------------------------------------------------------

ShiftedLat::ShiftedLat( const double &londeg, const double &latdeg ) : LonLat()
{
  setup(londeg,latdeg);
  set_typeinfo();
}

//------------------------------------------------------------------------------

void ShiftedLat::setup(const eckit::Parametrisation& p)
{
  long nlon, nlat;

  if( p.get("N",N_ ) ) // --> global grid (2*N x N)
  {
    setup(N_);
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
  }
}

//------------------------------------------------------------------------------

void ShiftedLat::setup( const long N )
{
  double delta = 90./static_cast<double>(N);
  std::vector<double> lats(2*N);
  std::vector<long>   nlons(2*N,4*N);
  std::vector<double> lonmin(2*N,0.);
  std::vector<double> lonmax(2*N,360.-delta);

  double latmax = 90.-0.5*delta;

  for( size_t jlat=0; jlat<2*N; ++jlat )
  {
    lats[jlat] = latmax - static_cast<double>(jlat)*delta;
  }

  ReducedGrid::N_ = N;
  ReducedGrid::setup(2*N,lats.data(),nlons.data(),lonmin.data(),lonmax.data());
}

//------------------------------------------------------------------------------

void ShiftedLat::setup(const long nlon, const long nlat)
{
  double londeg = 360./static_cast<double>(nlon);
  double latdeg = 180./static_cast<double>(nlat);

  std::vector<double> lats(nlat);
  std::vector<long>   nlons(nlat,nlon);
  std::vector<double> lonmin(nlat,0.);
  std::vector<double> lonmax(nlon,360.-londeg);

  double latmax = 90.-0.5*latdeg;

  for( size_t jlat=0; jlat<nlat; ++jlat )
  {
    lats[jlat] = latmax - static_cast<double>(jlat)*latdeg;
  }

  if( nlat%2 == 0 && nlon==2*nlat )
  {
    ReducedGrid::N_ = nlat/2;
  }
  ReducedGrid::setup(nlat,lats.data(),nlons.data(),lonmin.data(),lonmax.data());
}

void ShiftedLat::setup( const double londeg, const double latdeg )
{
  double Llon = 360.-londeg;
  double Llat = 180.-latdeg;
  double nlon_real = Llon/londeg + 1.;
  double nlat_real = Llat/latdeg + 1.;
  long nlon = static_cast<long>(nlon_real);
  long nlat = static_cast<long>(nlat_real);
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

eckit::Properties ShiftedLat::spec() const
{
  eckit::Properties grid_spec;

  grid_spec.set("grid_type",grid_type_str() );
  grid_spec.set("N", N() );

  BoundBox bbox = boundingBox();
  grid_spec.set("bbox_s", bbox.min().lat());
  grid_spec.set("bbox_w", bbox.min().lon());
  grid_spec.set("bbox_n", bbox.max().lat());
  grid_spec.set("bbox_e", bbox.max().lon());

  return grid_spec;
}

//-----------------------------------------------------------------------------

} // namespace global
} // namespace lonlat
} // namespace grid
} // namespace atlas
