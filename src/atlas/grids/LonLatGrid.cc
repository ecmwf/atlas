/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <typeinfo> // std::bad_cast

#include "eckit/memory/Builder.h"

#include "atlas/atlas_config.h"
#include "atlas/grids/LonLatGrid.h"

using eckit::BadParameter;
using eckit::Params;

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid,LonLatGrid,LonLatGrid::grid_type_str());

std::string LonLatGrid::grid_type_str()
{
  return "regular_ll";
}

std::string LonLatGrid::className()
{
  return "atlas.grid.LonLatGrid";
}

void LonLatGrid::set_typeinfo()
{
  std::stringstream s;
  s << "ll." << nlon() << "x" << nlat();
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

LonLatGrid::LonLatGrid() : ReducedLonLatGrid()
{
}

LonLatGrid::LonLatGrid(const eckit::Parametrisation& p)
{
  setup(p);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const size_t nlon, const size_t nlat, const BoundBox& bbox )
{
  setup(nlon,nlat,bbox);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const size_t nlon, const size_t nlat, TYPE poles )
{
  setup(nlon,nlat,poles);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const size_t nlat, TYPE poles )
{
  int nlon = 2*nlat;
  setup(nlon,nlat,poles);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const double londeg, const double latdeg, TYPE poles )
{
  setup(londeg,latdeg,poles);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const double londeg, const double latdeg, const BoundBox& bbox )
{
  setup(londeg,latdeg,bbox);
  set_typeinfo();
}


void LonLatGrid::setup(const eckit::Parametrisation& p)
{
  size_t nlon, nlat;

  bool poles(defaults::poles());
  p.get("poles",poles);

  if( p.get("N",N_ ) ) // --> global grid (2*N x N)
  {
    nlat = N_;
    if( poles )
      nlon = 2*(N_-1);
    else
      nlon = 2*N_;
    setup(nlon,nlat,poles);
  }
  else
  {
    if( !p.has("nlon") && !p.has("lon_inc") ) throw BadParameter("nlon or lon_inc missing in Params",Here());
    if( !p.has("nlat") && !p.has("lat_inc") ) throw BadParameter("nlat or lat_inc missing in Params",Here());

    bool bbox_present = p.has("bbox_n") && p.has("bbox_s") && p.has("bbox_e") && p.has("bbox_w");
    if( bbox_present ) // --> limited area grid
    {
      double bbox_n, bbox_s, bbox_e, bbox_w;
      p.get("bbox_n",bbox_n);
      p.get("bbox_s",bbox_s);
      p.get("bbox_e",bbox_e);
      p.get("bbox_w",bbox_w);
      BoundBox bbox(bbox_n,bbox_s,bbox_e,bbox_w);

      double lon_inc, lat_inc;
      if (p.get("nlon",nlon) && p.get("nlat",nlat))
      {
        setup(nlon,nlat,bbox);
      }
      else if (p.get("lon_inc",lon_inc) && p.get("lat_inc",lat_inc))
      {
        setup(lon_inc,lat_inc,bbox);
      }
      else
      {
        throw BadParameter("Bad combination of parameters");
      }
    }
    else // --> global grid (nlon x nlat)
    {
      double lon_inc, lat_inc;
      if (p.get("nlon",nlon) && p.get("nlat",nlat))
      {
        setup(nlon,nlat,poles);
      }
      else if (p.get("lon_inc",lon_inc) && p.get("lat_inc",lat_inc))
      {
        setup(lon_inc,lat_inc,poles);
      }
      else
      {
        throw BadParameter("Bad combination of parameters");
      }
    }
  }
}

void LonLatGrid::setup( const size_t nlon, const size_t nlat, const BoundBox& bbox )
{
  std::vector<double> lats(nlat);
  std::vector<long>   nlons(nlat,nlon);
  std::vector<double> lonmin(nlat,bbox.min().lon());
  std::vector<double> lonmax(nlat,bbox.max().lon());

  double latmin = bbox.min().lat();
  double latmax = bbox.max().lat();

  double delta = (latmax-latmin)/static_cast<double>(nlat-1);

  for( size_t jlat=0; jlat<nlat; ++jlat )
  {
    lats[jlat] = latmax - static_cast<double>(jlat)*delta;
  }

  ReducedGrid::setup(nlat,lats.data(),nlons.data(),lonmin.data(),lonmax.data());
}

void LonLatGrid::setup( const size_t nlon, const size_t nlat, bool poles )
{
    double londelta = 360./static_cast<double>(nlon);
    double latdelta = 180./static_cast<double>(nlat);

    BoundBox bbox = poles ?
                        BoundBox( 90.,-90, 360.-londelta, 0. ) :
                        BoundBox( 90.-0.5*latdelta, -90+0.5*latdelta, 360.-londelta, 0. );

    setup(nlon,nlat,bbox);
}


void LonLatGrid::setup( const double londeg, const double latdeg, bool poles )
{
  if( poles)
  {
    BoundBox bbox(90.,-90.,360.-londeg,0.);
    setup(londeg,latdeg,bbox);
  }
  else
  {
    BoundBox bbox(90.-0.5*latdeg,-90.+0.5*latdeg,360.-londeg,0.);
    setup(londeg,latdeg,bbox);
  }
}


void LonLatGrid::setup( const double londeg, const double latdeg, const BoundBox& bbox )
{
  double Llon = (bbox.max().lon()-bbox.min().lon());
  double Llat = (bbox.max().lat()-bbox.min().lat());
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
  setup(nlon,nlat,bbox);
}


eckit::Properties LonLatGrid::spec() const
{
  eckit::Properties grid_spec;

  grid_spec.set("grid_type",grid_type_str() );

  grid_spec.set("nlon", nlon() );
  grid_spec.set("nlat", nlat() );

  BoundBox bbox = boundingBox();
  grid_spec.set("bbox_s", bbox.min().lat());
  grid_spec.set("bbox_w", bbox.min().lon());
  grid_spec.set("bbox_n", bbox.max().lat());
  grid_spec.set("bbox_e", bbox.max().lon());

  return grid_spec;
}

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas
