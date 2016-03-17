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
#include "atlas/grid/global/gaussian/Gaussian.h"
#include "atlas/grid/global/gaussian/latitudes/Latitudes.h"

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid,ReducedGaussianGrid,ReducedGaussianGrid::grid_type_str());

std::string ReducedGaussianGrid::className()
{
  return "atlas.grid.ReducedGaussianGrid";
}

void ReducedGaussianGrid::set_typeinfo()
{
  std::stringstream s;
  s << "rgg.N" << N();
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

ReducedGaussianGrid::ReducedGaussianGrid()
{
}

ReducedGaussianGrid::ReducedGaussianGrid( const size_t N, const long nlons[], const Domain& d)
  : ReducedGrid(d)
{
  ReducedGrid::N_ = N;

  setup_N_hemisphere(N,nlons);
  set_typeinfo();
}

ReducedGaussianGrid::ReducedGaussianGrid(const eckit::Parametrisation& params)
{
  setup(params);
  set_typeinfo();
}

void ReducedGaussianGrid::setup( const eckit::Parametrisation& params )
{
  if( ! params.has("N") ) throw eckit::BadParameter("N missing in Params",Here());
  size_t N;
  params.get("N",N);
  ReducedGrid::N_ = N;

  if( ! params.has("latitudes") )
  {
    std::vector<long> nlons;
    params.get("npts_per_lat",nlons);
    setup_N_hemisphere(N,nlons.data());
  }
  else
  {
    ReducedGrid::setup(params);
  }
}

void ReducedGaussianGrid::setup_N_hemisphere( const size_t N, const long nlons[] )
{
  // hemisphere
  std::vector<double> lats (N);
  gaussian_latitudes_npole_equator(N,lats.data());
  ReducedGrid::setup_lat_hemisphere(N,lats.data(),nlons);
}

eckit::Properties ReducedGaussianGrid::spec() const
{
  eckit::Properties grid_spec;

  grid_spec.set("grid_type", grid_type_str() );

  grid_spec.set("nlat",nlat());
  grid_spec.set("N", N() );

  grid_spec.set("npts_per_lat",eckit::makeVectorValue(npts_per_lat()));
  grid_spec.set("latitudes",eckit::makeVectorValue(latitudes()));

  BoundBox bbox = boundingBox();
  grid_spec.set("bbox_s", bbox.min().lat());
  grid_spec.set("bbox_w", bbox.min().lon());
  grid_spec.set("bbox_n", bbox.max().lat());
  grid_spec.set("bbox_e", bbox.max().lon());

  return grid_spec;
}

//-----------------------------------------------------------------------------

extern "C" {

ReducedGrid* atlas__new_reduced_gaussian_grid(int nlon[], int nlat)
{
  std::vector<long> nlon_vector;
  nlon_vector.assign(nlon,nlon+nlat);
  return new ReducedGaussianGrid(nlat,nlon_vector.data());
}

}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace atlas
