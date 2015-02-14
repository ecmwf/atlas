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


#include "atlas/grids/ReducedGaussianGrid.h"
#include "atlas/grids/GaussianLatitudes.h"
#include "atlas/GridSpec.h"
#include "atlas/util/Debug.h"

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid,ReducedGaussianGrid,ReducedGaussianGrid::grid_type_str());

std::string ReducedGaussianGrid::className()
{
  return "atlas.grid.ReducedGaussianGrid";
}

void ReducedGaussianGrid::set_typeinfo()
{
  std::stringstream stream;
  stream << "rgg.N"<<N();
  uid_ = stream.str();
  grid_type_ = grid_type_str();
}

ReducedGaussianGrid::ReducedGaussianGrid() : ReducedGrid()
{
}

ReducedGaussianGrid::ReducedGaussianGrid( const int N, const int nlons[] )
{
  ReducedGrid::N_ = N;

  setup_N_hemisphere(N,nlons);
  set_typeinfo();
}

ReducedGaussianGrid::ReducedGaussianGrid(const eckit::Params& params)
{
  setup(params);
  mask(params);
  set_typeinfo();
}

void ReducedGaussianGrid::setup( const eckit::Params& params )
{
  if( ! params.has("N") ) throw eckit::BadParameter("N missing in Params",Here());
  int N = params.get("N");
  ReducedGrid::N_ = N;

  if( ! params.has("latitudes") )
  {
    eckit::ValueList list = params.get("npts_per_lat");
    std::vector<int> nlons(list.size());
    for(int j=0; j<nlons.size(); ++j)
      nlons[j] = list[j];

    setup_N_hemisphere(N,nlons.data());
  }
  else
  {
    ReducedGrid::setup(params);
  }
}

void ReducedGaussianGrid::setup_N_hemisphere( const int N, const int nlons[] )
{
  // hemisphere
  std::vector<double> lats (N);
  gaussian_latitudes_npole_equator(N,lats.data());
  ReducedGrid::setup_lat_hemisphere(N,lats.data(),nlons,DEG);
}



GridSpec ReducedGaussianGrid::spec() const
{
  GridSpec grid_spec( grid_type_str() );

  grid_spec.uid( uid() );
  grid_spec.set("nlat",nlat());
  grid_spec.set_npts_per_lat(npts_per_lat());
  grid_spec.set_latitudes(latitudes());
  grid_spec.set("N", N() );

  grid_spec.set_bounding_box(bounding_box());

  return grid_spec;
}

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas
