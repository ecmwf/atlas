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
#include <eckit/memory/Builder.h>

#include "atlas/atlas_config.h"
#ifdef ECKIT_HAVE_GRIB
#include "grib_api.h" // for grib_get_gaussian_latitudes()
#endif

#include "atlas/grids/ReducedGaussianGrid.h"
#include "atlas/GridSpec.h"
#include "atlas/util/Debug.h"

namespace atlas {
namespace grids {

std::string ReducedGaussianGrid::className()
{
  return "atlas.grid.ReducedGaussianGrid";
}

void ReducedGaussianGrid::set_typeinfo()
{
  std::stringstream stream;
  stream << "reduced_gg.N"<<N();
  uid_ = stream.str();
  hash_ = stream.str();
  grid_type_ = "reduced_gg";
}

ReducedGaussianGrid::ReducedGaussianGrid() : ReducedGrid()
{
}

ReducedGaussianGrid::ReducedGaussianGrid( const int N, const int nlons[] )
{
  setup_N_hemisphere(N,nlons);
  set_typeinfo();
}

ReducedGaussianGrid::ReducedGaussianGrid(const eckit::Params& params)
{
  setup(params);
  crop(params);
  set_typeinfo();
}

void ReducedGaussianGrid::setup( const eckit::Params& params )
{
  if( ! params.has("N") ) throw eckit::BadParameter("N missing in Params",Here());
  int N = params.get("N");

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
#ifdef ECKIT_HAVE_GRIB
  std::vector<double> lats (2*N);
  /// @todo this code should be moved into Atlas library and co-maintained with NA section
  grib_get_gaussian_latitudes(N, lats.data());
  ReducedGrid::setup_lat_hemisphere(N,nlons,lats.data(),DEG);
#else
  // hemisphere
  std::vector<double> lats (N);
  predict_gaussian_latitudes_hemisphere(N,lats.data());
  ReducedGrid::setup_lat_hemisphere(N,nlons,lats.data(),DEG);
#endif
}



GridSpec ReducedGaussianGrid::spec() const
{
  GridSpec grid_spec( ReducedGrid::spec() );

  grid_spec.set("N", N() );

  return grid_spec;
}

} // namespace grids
} // namespace atlas
