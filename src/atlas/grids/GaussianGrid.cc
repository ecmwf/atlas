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

#include "atlas/grids/GaussianGrid.h"
#include "atlas/grids/GaussianLatitudes.h"
#include "atlas/GridSpec.h"
#include "atlas/util/Debug.h"

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid, GaussianGrid,GaussianGrid::grid_type_str());

std::string GaussianGrid::className()
{
  return "atlas.grid.GaussianGrid";
}

void GaussianGrid::set_typeinfo()
{
  std::stringstream s;
  s << "gg.N"<< N();
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

GaussianGrid::GaussianGrid() : ReducedGaussianGrid()
{
}

GaussianGrid::GaussianGrid(const eckit::Params& params)
{
  eckit::ValueList list;

  if( ! params.has("N") ) throw eckit::BadParameter("N missing in Params",Here());

  int N = params["N"];

  N_ = N;

  if( ! params.has("latitudes") )
  {
    setup(N);
  }
  else
  {
    std::vector<int>    nlons(2*N,4*N);
    std::vector<double> lat;

    list = params["latitudes"];
    ASSERT(list.size() == 2*N);
    lat.resize( 2*N );
    for(int j=0; j<2*N; ++j)
      lat[j] = list[j];
    ReducedGrid::setup(lat.size(),lat.data(),nlons.data());
  }

  set_typeinfo();
}

GaussianGrid::GaussianGrid( const int N )
{
  ReducedGrid::N_ = N;
  setup(N);
  set_typeinfo();
}

void GaussianGrid::setup(const int N)
{
  std::vector<double> lats (N);
  gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data());
}

void GaussianGrid::setup_lat_hemisphere(const int N, const double lats[])
{
  std::vector<int> nlons(N,4*N);
  ReducedGrid::setup_lat_hemisphere(N,lats,nlons.data(),DEG);
}

GridSpec GaussianGrid::spec() const
{
  GridSpec grid_spec( grid_type_str() );

  grid_spec.set("N", N() );
  grid_spec.set_latitudes(latitudes());

  grid_spec.set("nlat",nlat());

  grid_spec.set_bounding_box(boundingBox());

  return grid_spec;
}

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas
