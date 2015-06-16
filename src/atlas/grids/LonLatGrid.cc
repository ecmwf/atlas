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

#include "eckit/memory/Builder.h"

#include "atlas/atlas_config.h"
#include "atlas/grids/LonLatGrid.h"
#include "atlas/GridSpec.h"
#include "atlas/util/Debug.h"

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

LonLatGrid::LonLatGrid(const Params& p)
{
  setup(p);
  mask(p);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const int nlon, const int nlat, const BoundBox& bbox )
{
  setup(nlon,nlat,bbox);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const int nlon, const int nlat, TYPE poles )
{
  setup(nlon,nlat,poles);
  set_typeinfo();
}

LonLatGrid::LonLatGrid( const int nlat, TYPE poles )
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


void LonLatGrid::setup(const Params& p)
{
  int nlon, nlat;

  bool poles(defaults::poles());
  if( p.has("poles") )
    poles = p["poles"];

  if( p.has("N") ) // --> global grid (2*N x N)
  {
    N_ = p["N"];
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
      BoundBox bbox(p["bbox_n"],p["bbox_s"],p["bbox_e"],p["bbox_w"]);
      if (p.has("nlon") && p.has("nlat"))
      {
        int nlon = p["nlon"];
        int nlat = p["nlat"];
        setup(nlon,nlat,bbox);
      }
      else if (p.has("lon_inc") && p.has("lat_inc"))
      {
        double londeg = p["lon_inc"];
        double latdeg = p["lat_inc"];
        setup(londeg,latdeg,bbox);
      }
    }
    else // --> global grid (nlon x nlat)
    {
      if (p.has("nlon") && p.has("nlat"))
      {
        int nlon = p["nlon"];
        int nlat = p["nlat"];
        setup(nlon,nlat,poles);
      }
      else if (p.has("lon_inc") && p.has("lat_inc"))
      {
        double londeg = p["lon_inc"];
        double latdeg = p["lat_inc"];
        setup(londeg,latdeg,poles);
      }
      else
      {
        throw BadParameter("Bad combination of parameters");
      }
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

void LonLatGrid::setup( const int nlon, const int nlat, bool poles )
{
  if( poles )
  {
    double londelta = 360./static_cast<double>(nlon);
    BoundBox bbox( 90.,-90, 360.-londelta, 0. );
    setup(nlon,nlat,bbox);
  }
  else
  {
    double latdelta = 180./static_cast<double>(nlat);
    double londelta = 360./static_cast<double>(nlon);
    BoundBox bbox( 90.-0.5*latdelta, -90+0.5*latdelta, 360.-londelta, 0. );
    setup(nlon,nlat,bbox);
  }
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
  int nlon = static_cast<int>(nlon_real);
  int nlat = static_cast<int>(nlat_real);
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


GridSpec LonLatGrid::spec() const
{
  GridSpec grid_spec( grid_type_str() );

  grid_spec.set("nlon", nlon() );
  grid_spec.set("nlat", nlat() );

  grid_spec.set_bounding_box(bounding_box());

  return grid_spec;
}

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas
