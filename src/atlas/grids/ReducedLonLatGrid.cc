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
#include "atlas/grids/ReducedLonLatGrid.h"
#include "atlas/GridSpec.h"
#include "atlas/util/Debug.h"

using eckit::Params;
using eckit::BadParameter;

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid,ReducedLonLatGrid,ReducedLonLatGrid::grid_type_str());

std::string ReducedLonLatGrid::grid_type_str()
{
  return "reduced_ll";
}

std::string ReducedLonLatGrid::className()
{
  return "atlas.grid.ReducedLonLatGrid";
}

void ReducedLonLatGrid::set_typeinfo()
{
  std::stringstream s;
  s << "rll.N" << N();
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

ReducedLonLatGrid::ReducedLonLatGrid()
  : ReducedGrid()
{
}

ReducedLonLatGrid::ReducedLonLatGrid( const int nlat, const int nlons[], bool poles, const Domain& domain)
 : ReducedGrid(domain)
{
  ReducedGrid::N_ = nlat;
  poles_ = poles;
  setup(nlat,nlons,poles_);
  set_typeinfo();
}

ReducedLonLatGrid::ReducedLonLatGrid( const Params& params )
{
  setup(params);
  set_typeinfo();
}

void ReducedLonLatGrid::setup( const Params& params )
{
  if( ! params.has("nlat") )         throw BadParameter("N missing in Params",Here());
  if( ! params.has("npts_per_lat") ) throw BadParameter("npts_per_lat missing in Params",Here());

  int nlat = params["nlat"];

  if( params.has("N") )
    N_ = params["N"];
  else
    N_ = nlat;

  poles_ = defaults::poles();
  if( params.has("poles") ) poles_ = params["poles"];


  eckit::ValueList list = params["npts_per_lat"];
  std::vector<int> nlons(list.size());
  for(int j=0; j<nlons.size(); ++j)
    nlons[j] = list[j];

  if( params.has("latitudes") )
  {
    ReducedGrid::setup(params);
  }
  else
  {
    setup(nlat,nlons.data(),poles_);
  }
}

void ReducedLonLatGrid::setup( const int nlat, const int nlons[], bool poles )
{
  std::vector<double> lats (nlat);

  double delta = domain_.north() - domain_.south();
  double latmax;

  if( poles )
  {
    delta = delta / static_cast<double>(nlat-1);
    latmax = domain_.north();
  }
  else
  {
    delta = delta / static_cast<double>(nlat);
    latmax = domain_.north() - 0.5*delta;
  }

  for( int jlat=0; jlat<nlat; ++jlat )
  {
    lats[jlat] = latmax - static_cast<double>(jlat)*delta;
  }
  ReducedGrid::setup(lats.size(),lats.data(),nlons);
}


GridSpec ReducedLonLatGrid::spec() const
{
  GridSpec grid_spec( grid_type_str() );

  grid_spec.set("N", N() );
  grid_spec.set("nlat", nlat() );
  grid_spec.set_npts_per_lat(npts_per_lat());

  if( nlat() != N() )
    grid_spec.set_latitudes(latitudes());

  grid_spec.set("poles",poles_);
  grid_spec.set_bounding_box(boundingBox());

  return grid_spec;
}

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas
