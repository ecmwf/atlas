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
#include <eckit/filesystem/PathName.h>
#include <eckit/memory/Factory.h>
#include <eckit/memory/Builder.h>
#include <eckit/parser/JSON.h>

#include "atlas/atlas.h"
#include "atlas/GridSpec.h"
#include "atlas/GridSpecParams.h"
#include "atlas/grids/grids.h"


//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;
using namespace atlas::grids;

//------------------------------------------------------------------------------------------------------

class AtlasGrids : public eckit::Tool {

  virtual void run();

public:

  AtlasGrids(int argc,char **argv): eckit::Tool(argc,argv)
  {
    //atlas_init(argc,argv);
    do_run = false;

    bool help = Resource< bool >("-h",false);

    std::string help_str =
        "NAME\n"
        "       atlas-grids - Catalogue of available built-in grids\n"
        "\n"
        "SYNOPSIS\n"
        "       atlas-grids GRID [OPTION]...\n"
        "\n"
        "DESCRIPTION\n"
        "       Browse catalogue of grids\n"
        "\n"
        "       -a       List all grids. The names are the GRID argument\n"
        "\n"
        "       -i       List information about GRID\n"
        "\n"
        "       -json    Export json\n"
        "\n"
        "       -rtable  Export IFS rtable\n"
        "\n"
        "AUTHOR\n"
        "       Written by Willem Deconinck.\n"
        "\n"
        "ECMWF                        November 2014"
        ;

    if( help )
    {
      std::cout << help_str << std::endl;
      return;
    }

    key = "";
    for( int i=0; i<argc; ++i )
    {
      if( i==1 && argv[i][0] != '-' )
      {
        key = std::string(argv[i]);
      }
    }

    info = Resource< bool >("-i",false);
    if( info && !key.empty() ) do_run = true;

    json = Resource< bool >("-json",false);
    if( json && !key.empty() ) do_run = true;

    rtable = Resource< bool >("-rtable",false);
    if( rtable && !key.empty() ) do_run = true;

    all = false;
    if( do_run == false && key.empty() )
    {
      all = true;
      do_run = true;
    }
    if( !key.empty() && do_run == false )
    {
      Log::error() << "Option wrong or missing after '" << key << "'" << std::endl;
    }
  }

private:

  bool all;
  std::string key;
  bool info;
  bool json;
  bool rtable;
  bool do_run;
};

//------------------------------------------------------------------------------------------------------

void AtlasGrids::run()
{
  if( !do_run ) return;

  if( all )
  {
    std::vector<std::string> keys = Factory<Grid>::instance().keys();
    Log::info() << "Available grids:" << std::endl;
    for( int i=0; i< keys.size(); ++i )
    {
      Log::info() << "  -- " << keys[i] << std::endl;
    }
  }

  if( !key.empty() )
  {
    if( !Factory<Grid>::instance().exists(key) )
    {
      Log::error() << "Grid " << key << " was not found." << std::endl;
      return;
    }

    ReducedGrid::Ptr grid ( ReducedGrid::create(key) );
    if( info )
    {
      Log::info() << "Grid " << key << std::endl;
      Log::info() << "   type:                             "
                  << grid->grid_type() << std::endl;
      Log::info() << "   uid:                              "
                  << grid->uid() << std::endl;
      if( grid->grid_type() == GaussianGrid::gtype() )
      {
        Log::info() << "   N number:                         "
                    << dynamic_cast<GaussianGrid*>(grid.get())->N() << std::endl;
      }
      if( grid->grid_type() == ReducedGaussianGrid::gtype() )
      {
        Log::info() << "   N number:                         "
                    << dynamic_cast<ReducedGaussianGrid*>(grid.get())->N() << std::endl;
      }
      Log::info() << "   number of points:                 "
                  << grid->npts() << std::endl;
      Log::info() << "   number of latitudes (N-S):        "
                  << grid->nlat() << std::endl;
      Log::info() << "   number of longitudes (max):       "
                  << grid->nlonmax() << std::endl;
      Log::info() << "   approximate resolution:           "
                  << 180./grid->nlat() << " deg" << std::endl;
      Log::info() << "   approximate resolution for Earth: "
                  << 40075 / (2.* grid->nlat()) << " km" << std::endl;
      Log::info() << "   spectral truncation -- linear:    "
                  << grid->nlat() - 1 << std::endl;
      Log::info() << "   spectral truncation -- quadratic: "
                  << static_cast<int>(std::floor(2./3.*grid->nlat()+0.5))-1 << std::endl;
      Log::info() << "   spectral truncation -- cubic:     "
                  << static_cast<int>(std::floor(0.5*grid->nlat()+0.5))-1 << std::endl;
      Log::info() << "   bounding box -- lat N-S (deg):    "
                  << grid->bounding_box().max().lat() << ", "
                  << grid->bounding_box().min().lat() << std::endl;
      Log::info() << "   bounding box -- lon W-E (deg):    "
                  << grid->bounding_box().min().lon() << ", "
                  << grid->bounding_box().max().lon() << std::endl;
    }
    if( json )
    {
      std::stringstream stream;
      JSON js(stream);
      js.precision(16);
      js << grid->spec();
      std::cout << stream.str() << std::endl;
    }

    if( rtable )
    {
      std::stringstream stream;
      stream << "&NAMRGRI\n";
      for( int jlat=0; jlat<grid->nlat(); ++jlat )
        stream << " NRGRI("<< std::setfill('0') << std::setw(5) << 1+jlat <<")="<< std::setfill(' ') << std::setw(5) << grid->nlon(jlat) <<",\n";
      stream << "/" << std::flush;
      std::cout << stream.str() << std::endl;
    }
  }
  //atlas_finalize();
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  AtlasGrids tool(argc,argv);
  tool.start();
  return 0;
}
