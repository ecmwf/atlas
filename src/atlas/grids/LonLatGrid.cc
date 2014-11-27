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

#include "atlas/atlas_config.h"
#ifdef ECKIT_HAVE_GRIB
#include "grib_api.h" // for grib_get_LonLat_latitudes()
#endif

#include "atlas/grids/LonLatGrid.h"
#include "atlas/GridSpec.h"
#include "atlas/util/Debug.h"

using eckit::BadParameter;
using eckit::Params;

namespace atlas {
namespace grids {

std::string LonLatGrid::gtype()
{
  return "regular_ll";
}

std::string LonLatGrid::className()
{
  return "atlas.grid.LonLatGrid";
}

void LonLatGrid::set_typeinfo()
{
  std::stringstream stream;
  stream << gtype()<<"."<<nlon()<<"x"<<nlat();
  uid_ = stream.str();
  hash_ = stream.str();
  grid_type_ = gtype();
}

LonLatGrid::LonLatGrid() : ReducedGrid()
{
}

LonLatGrid::LonLatGrid(const Params& p)
{
  setup(p);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const int nlon, const int nlat, const BoundBox& bbox )
{
  setup(nlon,nlat,bbox);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const int nlon, const int nlat )
{
  setup(nlon,nlat);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const int N )
{
  int nlon = 4*N;
  int nlat = 2*N;
  setup(nlon,nlat);
  set_typeinfo();
}

void LonLatGrid::setup(const Params& p)
{
  int nlon, nlat;

  if( p.has("N") ) // --> global grid (4*N x 2*N)
  {
    const int N = p["N"];
    nlon = 4*N;
    nlat = 2*N;
    setup(nlon,nlat);
  }
  else
  {
    if( ! p.has("nlon") ) throw BadParameter("nlon missing in Params",Here());
    if( ! p.has("nlat") ) throw BadParameter("nlat missing in Params",Here());
    nlon = p["nlon"];
    nlat = p["nlat"];

    if( p.has("bbox_n") && p.has("bbox_s") && p.has("bbox_e") && p.has("bbox_w") )
    {
      // --> limited area grid
      setup(nlon,nlat,BoundBox(p["bbox_n"],p["bbox_s"],p["bbox_e"],p["bbox_w"]));
    }
    else // --> global grid (nlon x nlat)
    {
      setup(nlon,nlat);
    }
  }
}

void LonLatGrid::setup( const int nlon, const int nlat, const BoundBox& bbox )
{
  std::vector<double> lats(nlat);
  std::vector<int>    nlons(nlat,nlon);
  std::vector<double> lonmin(nlat,bbox.min().lon());
  std::vector<double> lonmax(nlat,bbox.max().lon());

  double latmin = bbox.min().lat();
  double latmax = bbox.max().lat();

  double delta = (latmax-latmin)/static_cast<double>(nlat-1);

  for( int jlat=0; jlat<nlat; ++jlat )
  {
    lats[jlat] = latmax - static_cast<double>(jlat)*delta;
  }

  ReducedGrid::setup(nlat,lats.data(),nlons.data(),lonmin.data(),lonmax.data());
}

void LonLatGrid::setup(const int nlon, const int nlat)
{
  double latdelta = 180./static_cast<double>(nlat);
  double londelta = 360./static_cast<double>(nlon);
  BoundBox bbox( 90.-0.5*latdelta, -90+0.5*latdelta, 360.-londelta, 0. );
  setup(nlon,nlat,bbox);
}

GridSpec LonLatGrid::spec() const
{
  GridSpec grid_spec( gtype() );

  grid_spec.set("hash", hash() );
  grid_spec.set("uid", uid()  );

  grid_spec.set("nlon", nlon() );
  grid_spec.set("nlat", nlat() );

  grid_spec.set_bounding_box(bounding_box());

  return grid_spec;
}

} // namespace grids
} // namespace atlas
