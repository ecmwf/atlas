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

#include "atlas/grids/GaussianGrid.h"
#include "atlas/GridSpec.h"
#include "atlas/util/Debug.h"

namespace atlas {
namespace grids {

std::string GaussianGrid::className()
{
  return "atlas.grid.GaussianGrid";
}

void GaussianGrid::set_typeinfo()
{
  std::stringstream stream;
  stream << "regular_gg.N"<<N();
  uid_ = stream.str();
  hash_ = stream.str();
  grid_type_ = "regular_gg";
}

GaussianGrid::GaussianGrid() : ReducedGrid()
{
}

GaussianGrid::GaussianGrid(const eckit::Params& params)
{
  eckit::ValueList list;

  if( ! params.has("N") ) throw eckit::BadParameter("N missing in Params",Here());

  int N = params.get("N");

  if( ! params.has("latitudes") )
  {
    setup(N);
  }
  else
  {
    std::vector<int>    nlons(2*N,4*N);
    std::vector<double> lat;

    list = params.get("latitudes");
    ASSERT(list.size() == 2*N);
    lat.resize( 2*N );
    for(int j=0; j<2*N; ++j)
      lat[j] = list[j];
    ReducedGrid::setup(lat.size(),lat.data(),nlons.data());
  }

  crop(params);
  set_typeinfo();
}

GaussianGrid::GaussianGrid( const int N )
{
  setup(N);
  set_typeinfo();
}

void GaussianGrid::setup(const int N)
{
#ifdef ECKIT_HAVE_GRIB
  std::vector<int>    nlons(2*N,4*N);
  std::vector<double> lats (2*N);

  /// @todo this code should be moved into Atlas library and co-maintained with NA section
  grib_get_gaussian_latitudes(N, lats.data());
  ReducedGrid::setup(2*N,lats.data(),nlons.data());
#else
  // hemisphere
  std::vector<double> lats (N);
  predict_gaussian_latitudes_hemisphere(N,lats.data());
  setup_lat_hemisphere(N,lats.data());
#endif
}

void GaussianGrid::setup_lat_hemisphere(const int N, const double lats[])
{
  std::vector<int> nlons(N,4*N);
  ReducedGrid::setup_lat_hemisphere(N,lats,nlons.data(),DEG);
}

GridSpec GaussianGrid::spec() const
{
  GridSpec grid_spec( ReducedGrid::spec() );

  grid_spec.set("N", N() );

  return grid_spec;
}

} // namespace grids
} // namespace atlas
