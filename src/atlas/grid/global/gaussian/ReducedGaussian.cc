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
#include "atlas/grid/global/gaussian/ReducedGaussian.h"

namespace atlas {
namespace grid {
namespace global {
namespace gaussian {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid,ReducedGaussian,ReducedGaussian::grid_type_str());

std::string ReducedGaussian::className()
{
  return "atlas.grid.global.gaussian.ReducedGaussian";
}

void ReducedGaussian::set_typeinfo()
{
  std::stringstream s;
  s << "reduced_gaussian.N" << N();
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

ReducedGaussian::ReducedGaussian( const size_t N, const long nlons[] )
  : Gaussian()
{
  Structured::N_ = N;

  setup_N_hemisphere(N,nlons);
  set_typeinfo();
}

ReducedGaussian::ReducedGaussian(const eckit::Parametrisation& params)
  : Gaussian()
{
  setup(params);
  set_typeinfo();
}

void ReducedGaussian::setup( const eckit::Parametrisation& params )
{
  if( ! params.has("N") ) throw eckit::BadParameter("N missing in Params",Here());
  size_t N;
  params.get("N",N);
  Structured::N_ = N;

  if( ! params.has("latitudes") )
  {
    std::vector<long> nlons;
    params.get("npts_per_lat",nlons);
    setup_N_hemisphere(N,nlons.data());
  }
  else
  {
    Structured::setup(params);
  }
}

eckit::Properties ReducedGaussian::spec() const
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

Structured* atlas__new_reduced_gaussian_grid(int nlon[], int nlat)
{
  std::vector<long> nlon_vector;
  nlon_vector.assign(nlon,nlon+nlat);
  return new ReducedGaussian(nlat,nlon_vector.data());
}

}

//-----------------------------------------------------------------------------

} // namespace gaussian
} // namespace global
} // namespace grid
} // namespace atlas
