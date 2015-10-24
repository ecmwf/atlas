/*
 * (C) Copyright 1996-2015 ECMWF.
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
  s << "F"<< N();
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

GaussianGrid::GaussianGrid() : ReducedGaussianGrid()
{
}

GaussianGrid::GaussianGrid(const eckit::Parametrisation& params)
{
  if( ! params.get("N",N_) ) throw eckit::BadParameter("N missing in Params",Here());

  if( ! params.has("latitudes") )
  {
    setup(N_);
  }
  else
  {
    std::vector<long>    nlons(2*N_,4*N_);
    std::vector<double> lat;

    params.get("latitudes",lat);
    ASSERT(lat.size() == 2*N_);
    ReducedGrid::setup(lat.size(),lat.data(),nlons.data());
  }

  set_typeinfo();
}

GaussianGrid::GaussianGrid( const size_t N )
{
  ReducedGrid::N_ = N;
  setup(N);
  set_typeinfo();
}

void GaussianGrid::setup(const size_t N)
{
  std::vector<double> lats (N);
  gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data());
}

void GaussianGrid::setup_lat_hemisphere(const size_t N, const double lats[])
{
  std::vector<long> nlons(N,4*N);
  ReducedGrid::setup_lat_hemisphere(N,lats,nlons.data(),DEG);
}

eckit::Properties GaussianGrid::spec() const
{
  eckit::Properties grid_spec;
  grid_spec.set("grid_type",grid_type_str());

  grid_spec.set("N", N() );
  grid_spec.set("latitudes",eckit::makeVectorValue(latitudes()));

  grid_spec.set("nlat",nlat());

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
