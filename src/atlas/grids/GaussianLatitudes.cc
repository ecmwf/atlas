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

#include "atlas/grids/GaussianLatitudes.h"

#include "eckit/memory/ScopedPtr.h"

#include "atlas/atlas_config.h"
#include "atlas/Parameters.h"
#include "atlas/grids/gausslat/gausslat.h"

#ifdef ECKIT_HAVE_GRIB
  #include "grib_api.h"
#endif

using eckit::ConcreteBuilderT0;
using eckit::Factory;
using eckit::ScopedPtr;
using eckit::Log;

using atlas::grids::gausslat::GaussianLatitudes;

namespace atlas {
namespace grids {

//----------------------------------------------------------------------------------------------------------------------

void predict_gaussian_colatitudes_hemisphere(const size_t N, double colat[]);

void predict_gaussian_latitudes_hemisphere(const size_t N, double lat[]);

namespace {
  
void colat_to_lat_hemisphere(const size_t N, const double colat[], double lats[], const AngleUnit unit)
{
  std::copy( colat, colat+N, lats );
  double pole = (unit == DEG ? 90. : M_PI_2);
  for(size_t i=0; i<N; ++i)
    lats[i]=pole-lats[i];
}

} //  anonymous namespace

//----------------------------------------------------------------------------------------------------------------------

void gaussian_latitudes_npole_equator(const size_t N, double lats[])
{
  std::stringstream Nstream; Nstream << N;
  std::string Nstr = Nstream.str();
  if( Factory<GaussianLatitudes>::instance().exists(Nstr) )
  {
    ScopedPtr<GaussianLatitudes> gl ( Factory<GaussianLatitudes>::instance().get(Nstr).create() );
    gl->assign(lats);
  }
  else
  {
#ifndef ECKIT_HAVE_GRIB
    Log::warning() << "Unfortunately, computation of the gaussian latitudes depends at the moment on grib_api.\n"
                   << "Atlas was built without grib_api support (via eckit)\n"
                   << "Resorting to use predicted latitudes instead (accuracy of 2 decimals)" << std::endl;

    predict_gaussian_latitudes_hemisphere(N,lats);
#else
    Log::warning() << "Using grib_api to compute gaussian latitudes..." << std::endl;
    std::vector<double> lats2(2*N);
    grib_get_gaussian_latitudes(N, lats2.data());
    for(size_t jlat=0; jlat<N; ++jlat )
      lats[jlat] = lats2[jlat];
#endif
  }
}

void gaussian_latitudes_npole_spole(const size_t N, double lats[])
{
    gaussian_latitudes_npole_equator(N,lats);
    size_t end = 2*N-1;
    for(size_t jlat=0; jlat<N; ++jlat) {
        lats[end--] = -lats[jlat];
    }
}

void predict_gaussian_colatitudes_hemisphere(const size_t N, double colat[])
{
  double z;
  for(size_t i=0; i<N; ++i )
  {
    z = (4.*(i+1.)-1.)*M_PI/(4.*2.*N+2.);
    colat[i] = ( z+1./(tan(z)*(8.*(2.*N)*(2.*N))) ) * Constants::radianToDegrees();
  }
}

void predict_gaussian_latitudes_hemisphere(const size_t N, double lats[])
{
  std::vector<double> colat(N);
  predict_gaussian_colatitudes_hemisphere(N,colat.data());
  colat_to_lat_hemisphere(N,colat.data(),lats,DEG);
}

//----------------------------------------------------------------------------------------------------------------------

// TODO: evaluate this, and possibly use
namespace {

  void initialGaussianLatitudes(size_t resolution, std::vector<double>& lats);
  void refineGaussianLatitude(size_t resolution, double& value);
  void gaussian(size_t resolution, std::vector<double>& lats);


  void gaussian(size_t resolution, std::vector<double>& lats)
  {
      lats.resize(resolution);

      // get first guess
      initialGaussianLatitudes(resolution, lats);

      // refine these values
      for (size_t i = 0; i < lats.size(); ++i)
          refineGaussianLatitude(resolution, lats[i]);
  }

  void refineGaussianLatitude(size_t resolution, double& value)
  {
      // NB code only tested on positive (Northern) hemisphere
      //
      static const double rad2deg = Constants::radianToDegrees();
      static const double convval = 1.0 - ((2.0 / M_PI)*(2.0 / M_PI)) * 0.25;
      static const double convergeance_precision = 1.0e-14;

      const size_t nlat = 2 * resolution;

      double root = cos(value / sqrt( ((((double)(nlat))+0.5)*(((double)(nlat))+0.5)) + convval) );

      double conv = 1.0;
      double mem1, mem2, legfonc;

      // Newton Iteration
      while(fabs(conv) > convergeance_precision)
      {
          mem2 = 1.0;
          mem1 = root;

          //  Compute Legendre polynomial
          for(size_t legi = 0; legi < nlat; legi++)
          {
              legfonc = ( (2.0 * (legi +1.0) - 1.0) * root * mem1 - legi * mem2) / ((double)(legi+1.0));
              mem2 = mem1;
              mem1 = legfonc;
          }

          //  Perform Newton iteration
          conv = legfonc / ((((double)(nlat)) * (mem2 - root * legfonc) ) / (1.0 - (root *root)));
          root -= conv;

          //  Routine fails if no convergence after JPMAXITER iterations.
      }

      //   Set North and South values using symmetry.

      value = asin(root) * rad2deg;
  }


  void initialGaussianLatitudes(size_t resolution, std::vector<double>& lats)
  {

      // Computes initial approximations for Gaussian latitudes
      // return zeroes of the Bessel function

      static const double zeroes[] =  {     2.4048255577E0,   5.5200781103E0,
                  8.6537279129E0,   11.7915344391E0,  14.9309177086E0,
                  18.0710639679E0,  21.2116366299E0,  24.3524715308E0,
                  27.4934791320E0,  30.6346064684E0,  33.7758202136E0,
                  36.9170983537E0,  40.0584257646E0,  43.1997917132E0,
                  46.3411883717E0,  49.4826098974E0,  52.6240518411E0,
                  55.7655107550E0,  58.9069839261E0,  62.0484691902E0,
                  65.1899648002E0,  68.3314693299E0,  71.4729816036E0,
                  74.6145006437E0,  77.7560256304E0,  80.8975558711E0,
                  84.0390907769E0,  87.1806298436E0,  90.3221726372E0,
                  93.4637187819E0,  96.6052679510E0,  99.7468198587E0,
                  102.8883742542E0, 106.0299309165E0, 109.1714896498E0,
                  112.3130502805E0, 115.4546126537E0, 118.5961766309E0,
                  121.7377420880E0, 124.8793089132E0, 128.0208770059E0,
                  131.1624462752E0, 134.3040166383E0, 137.4455880203E0,
                  140.5871603528E0, 143.7287335737E0, 146.8703076258E0,
                  150.0118824570E0, 153.1534580192E0, 156.2950342685E0 };

      ASSERT(sizeof(zeroes) > 0);

      const size_t nz = sizeof(zeroes) / sizeof(zeroes[0]);

      lats.resize(resolution);

      for(size_t i = 0; i < lats.size(); i++)
      {
          if(i < nz)
              lats[i] = zeroes[i];
          else
          {
              ASSERT(i > 0);
              lats[i] = lats[i-1] + M_PI;
          }
      }
  }

} // anonymous namespace

//----------------------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

