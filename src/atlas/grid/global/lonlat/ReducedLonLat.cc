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
#include "atlas/grid/global/lonlat/ReducedLonLat.h"

using namespace eckit;

namespace atlas {
namespace grid {
namespace global {
namespace lonlat {

//------------------------------------------------------------------------------

register_BuilderT1(Grid,ReducedLonLat,ReducedLonLat::grid_type_str());

std::string ReducedLonLat::grid_type_str()
{
  return "reduced_lonlat";
}

std::string ReducedLonLat::className()
{
  return "atlas.grid.global.lonlat.ReducedLonLat";
}

void ReducedLonLat::set_typeinfo()
{
  std::stringstream s;
  s << "reduced_lonlat";
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

ReducedLonLat::ReducedLonLat( 
    const size_t nlat,
    const long nlon[],
    const Shift& shift )
 : LonLat(Shift::NONE)
{
  setup(nlat,nlon,shift);
  set_typeinfo();
}

ReducedLonLat::ReducedLonLat( const eckit::Parametrisation& params )
  : LonLat(Shift::NONE)
{
  setup(params);
  set_typeinfo();
}

void ReducedLonLat::setup( const eckit::Parametrisation& params )
{
  if( ! params.has("nlat") )         throw BadParameter("N missing in Params",Here());
  if( ! params.has("npts_per_lat") ) throw BadParameter("npts_per_lat missing in Params",Here());

  size_t nlat;
  params.get("nlat",nlat);
  
  bool shift_lon = false;
  bool shift_lat = false;
  params.get("shift_lon",shift_lon);
  params.get("shift_lat",shift_lat);

  shift_ = Shift(shift_lon,shift_lat);

  std::vector<long> nlon;
  params.get("npts_per_lat",nlon);
  setup(nlat,nlon.data(),shift_);

  params.get("N",N_);
}

void ReducedLonLat::setup( 
    const size_t nlat, 
    const long nlon[],
    const Shift& shift )
{
  shift_ = shift;
  std::vector<double> lats (nlat);

  double delta = 180.;
  double latmax;

  if( shifted().lat() == false )
  {
    delta = delta / static_cast<double>(nlat-1);
    latmax = 90.;
  }
  else
  {
    delta = delta / static_cast<double>(nlat);
    latmax = 90. - 0.5*delta;
  }

  for(size_t jlat = 0; jlat < nlat; ++jlat)
  {
    lats[jlat] = latmax - static_cast<double>(jlat)*delta;
  }
  Structured::setup(lats.size(),lats.data(),nlon);
}


eckit::Properties ReducedLonLat::spec() const
{
  eckit::Properties grid_spec;

  grid_spec.set("grid_type", grid_type_str() );

  grid_spec.set("N", N() );
  grid_spec.set("nlat", nlat() );

  grid_spec.set("npts_per_lat",eckit::makeVectorValue(npts_per_lat()));

  grid_spec.set("latitudes",eckit::makeVectorValue(latitudes()));

  BoundBox bbox = boundingBox();
  grid_spec.set("bbox_s", bbox.min().lat());
  grid_spec.set("bbox_w", bbox.min().lon());
  grid_spec.set("bbox_n", bbox.max().lat());
  grid_spec.set("bbox_e", bbox.max().lon());
  grid_spec.set("shift_lon",shifted().lon());
  grid_spec.set("shift_lat",shifted().lat());

  return grid_spec;
}

//------------------------------------------------------------------------------

} // namespace lonlat
} // namespace global
} // namespace grid
} // namespace atlas
