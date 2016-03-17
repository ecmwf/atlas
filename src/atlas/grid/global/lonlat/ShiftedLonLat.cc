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
#include "atlas/grid/global/lonlat/ShiftedLonLat.h"

using eckit::BadParameter;
using eckit::Params;

namespace atlas {
namespace grid {
namespace global {
namespace lonlat {

//------------------------------------------------------------------------------

register_BuilderT1(Grid,ShiftedLonLat,ShiftedLonLat::grid_type_str());

//------------------------------------------------------------------------------

std::string ShiftedLonLat::grid_type_str()
{
  return "shifted_lonlat";
}

//------------------------------------------------------------------------------

std::string ShiftedLonLat::className()
{
  return "atlas.grid.global.lonlat.ShiftedLonLat";
}

//------------------------------------------------------------------------------

void ShiftedLonLat::set_typeinfo()
{
  std::stringstream s;
  s << "S" << N();
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

//------------------------------------------------------------------------------

ShiftedLonLat::ShiftedLonLat() : ReducedLonLatGrid()
{
}

//------------------------------------------------------------------------------

ShiftedLonLat::ShiftedLonLat(const eckit::Parametrisation& p)
{
  setup(p);
  set_typeinfo();
}

//------------------------------------------------------------------------------

ShiftedLonLat::ShiftedLonLat( const size_t N )
{
  setup(N);
  set_typeinfo();
}

//------------------------------------------------------------------------------

void ShiftedLonLat::setup(const eckit::Parametrisation& p)
{
  size_t nlon, nlat;

  if( p.get("N",N_ ) ) // --> global grid (2*N x N)
  {
    setup(N_);
  }
  else
  {
    throw BadParameter("\"N\" value missing.",Here());
  }
}

//------------------------------------------------------------------------------

void ShiftedLonLat::setup( const size_t N )
{
  double delta = 90./static_cast<double>(N);
  std::vector<double> lats(2*N);
  std::vector<long>   nlons(2*N,4*N);
  std::vector<double> lonmin(2*N,0.5*delta);
  std::vector<double> lonmax(2*N,360.-0.5*delta);

  double latmax = 90.-0.5*delta;

  for( size_t jlat=0; jlat<2*N; ++jlat )
  {
    lats[jlat] = latmax - static_cast<double>(jlat)*delta;
  }

  ReducedGrid::setup(2*N,lats.data(),nlons.data(),lonmin.data(),lonmax.data());
}

//------------------------------------------------------------------------------

eckit::Properties ShiftedLonLat::spec() const
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

//------------------------------------------------------------------------------

} // namespace global
} // namespace lonlat
} // namespace grid
} // namespace atlas
