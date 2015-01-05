/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <limits>
#include <cassert>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>

#include <eckit/exception/Exceptions.h>
#include <eckit/config/Resource.h>
#include <eckit/runtime/Tool.h>
#include <eckit/runtime/Context.h>

#include "atlas/atlas.h"
#ifdef ECKIT_HAVE_GRIB
#include "grib_api.h" // for grib_get_gaussian_latitudes()
#endif

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;

//------------------------------------------------------------------------------------------------------

class AtlasGaussianLatitudes : public eckit::Tool {

  virtual void run();

public:

  AtlasGaussianLatitudes(int argc,char **argv): eckit::Tool(argc,argv)
  {
    //atlas_init(argc,argv);
    do_run = false;

    bool help = Resource< bool >("-h",false);

    std::string help_str =
        "NAME\n"
        "       atlas-gaussian-latitudes - Compute gaussian latitudes, given the N number\n"
        "\n"
        "SYNOPSIS\n"
        "       atlas-gaussian-latitudes [OPTION]...\n"
        "\n"
        "DESCRIPTION\n"
        "       Compute gaussian latitudes, given the N number.\n"
        "       Latitudes start at the pole (+90), and decrease in value.\n"
        "\n"
        "       -N         Number of points between pole and equator\n"
        "\n"
        "       -full      If set, all latitudes will be given, otherwise only\n"
        "                  between North pole and equator.\n"
        "\n"
        "       -format    \"table\" (default)\n"
        "                  \"C\"\n"
        "\n"
        "       -compact   Write 5 latitudes per line if the format supports it\n"
        "\n"
        "AUTHOR\n"
        "       Written by Willem Deconinck.\n"
        "\n"
        "ECMWF                        December 2014"
        ;

    if( help )
    {
      std::cout << help_str << std::endl;
      return;
    }

    N = Resource< int >("-N",0);

    full = Resource< bool >("-full",false);

    compact = Resource< bool >("-compact",false);

    format = Resource< std::string >("-format", std::string("table") );

    if( N > 0 ) do_run = true;

#ifndef ECKIT_HAVE_GRIB
    std::cout << "Unfortunately, the computations of the gaussian latitudes require for now grib_api dependency" << std::endl;
    do_run = false;
#endif

  }

private:

  int N;
  bool full;
  bool compact;
  std::string format;
  bool do_run;
};

//------------------------------------------------------------------------------------------------------

void AtlasGaussianLatitudes::run()
{
  if( !do_run ) return;

  std::vector<double> lats (2*N);

#ifdef ECKIT_HAVE_GRIB
  grib_get_gaussian_latitudes(N, lats.data());
#endif

  int end = full ? 2*N : N;

  if( format == "table" )
  {
    for( int jlat=0; jlat<end; ++jlat )
      std::cout << std::setw(4) << jlat+1 << std::setw(17) << std::setprecision(12) << std::fixed << lats[jlat] << std::endl;
  }
  if( format == "C" )
  {
    std::cout << "double lat[] = {" << std::endl;
    for( int jlat=0; jlat<end; ++jlat )
    {
      std::cout << std::setw(16) << std::setprecision(12) << std::fixed << lats[jlat];
      if( jlat != end-1 ) std::cout << ",";
      if( (compact && (jlat+1)%5==0) || !compact || jlat == end-1 )
        std::cout << std::endl;
    }
    std::cout << "};" << std::endl;
  }
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  AtlasGaussianLatitudes tool(argc,argv);
  tool.start();
  return 0;
}
