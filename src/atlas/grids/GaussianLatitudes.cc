/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Jan 2014

#include "atlas/atlas_config.h"
#include "atlas/grids/gausslat/gausslat.h"

#ifdef ECKIT_HAVE_GRIB
  #include "grib_api.h"
#endif

using eckit::ConcreteBuilderT0;
using eckit::Factory;
using eckit::SharedPtr;
using eckit::Log;

using atlas::grids::gausslat::GaussianLatitudes;

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

void gaussian_latitudes_npole_equator(const int N, double lats[])
{
  std::stringstream Nstream; Nstream << N;
  std::string Nstr = Nstream.str();
  if( Factory<GaussianLatitudes>::instance().exists(Nstr) )
  {
    SharedPtr<GaussianLatitudes> gl ( Factory<GaussianLatitudes>::instance().get(Nstr).create() );
    gl->assign(lats);
  }
  else
  {
#ifndef ECKIT_HAVE_GRIB
    throw eckit::NotImplemented("Unfortunately, the computations of the gaussian latitudes require for now grib_api dependency, and no hard-coded version is provided",Here());

#else
    Log::info() << "using grib..." << std::endl;
    std::vector<double> lats2(2*N);
    grib_get_gaussian_latitudes(N, lats2.data());
    for( int jlat=0; jlat<N; ++jlat )
      lats[jlat] = lats2[jlat];
#endif
  }
}

void gaussian_latitudes_npole_spole(const int N, double lats[])
{
  gaussian_latitudes_npole_equator(N,lats);
  int end = 2*N-1;
  for( int jlat=0; jlat<N; ++jlat )
    lats[end-jlat] = lats[jlat];
}

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

