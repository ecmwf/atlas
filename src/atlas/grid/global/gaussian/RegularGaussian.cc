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
#include "atlas/grid/global/gaussian/RegularGaussian.h"

namespace atlas {
namespace grid {
namespace global {
namespace gaussian {

//-----------------------------------------------------------------------------

register_BuilderT1(Grid, RegularGaussian,RegularGaussian::grid_type_str());

std::string RegularGaussian::className()
{
  return "atlas.grid.global.gaussian.RegularGaussian";
}

void RegularGaussian::set_typeinfo()
{
  std::stringstream s;
  s << "F"<< N();
  shortName_ = s.str();
  grid_type_ = grid_type_str();
}

RegularGaussian::RegularGaussian(const eckit::Parametrisation& params) :
  Gaussian()
{
  if( ! params.get("N",N_) ) throw eckit::BadParameter("N missing in Params",Here());

  setup(N_);
  set_typeinfo();
}

RegularGaussian::RegularGaussian( const size_t N ) :
  Gaussian()
{
  Structured::N_ = N;
  setup(N);
  set_typeinfo();
}

void RegularGaussian::setup(const size_t N)
{
  std::vector<long> nlons(N,4*N);
  setup_N_hemisphere(N,nlons.data());
}

eckit::Properties RegularGaussian::spec() const
{
  eckit::Properties grid_spec;
  grid_spec.set("grid_type",grid_type_str());

  grid_spec.set("N", N() );

  grid_spec.set("nlat",nlat());

  BoundBox bbox = boundingBox();
  grid_spec.set("bbox_s", bbox.min().lat());
  grid_spec.set("bbox_w", bbox.min().lon());
  grid_spec.set("bbox_n", bbox.max().lat());
  grid_spec.set("bbox_e", bbox.max().lon());

  return grid_spec;
}

//-----------------------------------------------------------------------------

extern "C" {
  
Structured* atlas__grid__global__gaussian__RegularGaussian(size_t N)
{
  return new RegularGaussian(N);
}

}

//-----------------------------------------------------------------------------

} // namespace gaussian
} // namespace global
} // namespace grid
} // namespace atlas
