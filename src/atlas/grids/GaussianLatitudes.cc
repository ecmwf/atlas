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
#include "atlas/Util.h"

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

void predict_gaussian_colatitudes_hemisphere(const int N, double colat[]);

void predict_gaussian_latitudes_hemisphere(const int N, double lat[]);

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
    Log::warning() << "Unfortunately, the computations of the gaussian latitudes require for now grib_api dependency,\nand no hard-coded version is provided. Using predicted latitudes instead (accuracy of 2 decimals)" << std::endl;
    predict_gaussian_latitudes_hemisphere(N,lats);
#else
    Log::warning() << "Using grib_api to compute gaussian latitudes..." << std::endl;
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
  for( int jlat=0; jlat<N; ++jlat ) {
    lats[end--] = -lats[jlat];
  }
}


namespace {
static const double rad_to_deg = 180.*M_1_PI;
}

void predict_gaussian_colatitudes_hemisphere(const int N, double colat[])
{
  double z;
  for( int i=0; i<N; ++i )
  {
    z = (4.*(i+1.)-1.)*M_PI/(4.*2.*N+2.);
    colat[i] = ( z+1./(tan(z)*(8.*(2.*N)*(2.*N))) ) * rad_to_deg;
  }
}

void predict_gaussian_latitudes_hemisphere(const int N, double lats[])
{
  std::vector<double> colat(N);
  predict_gaussian_colatitudes_hemisphere(N,colat.data());
  colat_to_lat_hemisphere(N,colat.data(),lats,DEG);
}
//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

