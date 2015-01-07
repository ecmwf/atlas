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

#include "eckit/exception/Exceptions.h"
#include "eckit/config/Resource.h"
#include "eckit/runtime/Tool.h"
#include "eckit/runtime/Context.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/memory/Factory.h"
#include "eckit/memory/Builder.h"
#include "eckit/parser/JSON.h"
#include "eckit/parser/Tokenizer.h"

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

    bool help = Resource< bool >("--help",false);

    std::string help_str =
        "NAME\n"
        "       atlas-grids - Catalogue of available built-in grids\n"
        "\n"
        "SYNOPSIS\n"
        "       atlas-grids GRID [OPTION]... [--help] \n"
        "\n"
        "DESCRIPTION\n"
        "       Browse catalogue of grids\n"
        "\n"
        "       GRID: unique identifier for grid \n"
        "           Example values: rgg.N80, rgg.TL159, gg.N40, ll.128x64\n"
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

ReducedGrid::Ptr create_grid(const std::string& key)
{
  int N=0;
  int nlon=0;
  int nlat=0;
  Tokenizer tokenize(".");
  std::vector<std::string> tokens;
  tokenize(key,tokens);
  Translator<std::string,int> to_int;
  std::string uid = tokens[0];
  if( uid == "ll" ) uid = "regular_ll";
  if( uid == "gg" ) uid = "regular_gg";
  if( uid == "reduced_gg") uid = "rgg";

  if( tokens.size() > 1 )
  {
    if( tokens[1][0] == 'N' )
    {
      std::string Nstr(tokens[1],1,tokens[1].size()-1);
      N = to_int(Nstr);
    }
    else
    {
      std::vector<std::string> lonlat;
      Tokenizer tokenize_lonlat("x");
      tokenize_lonlat(tokens[1],lonlat);
      if( lonlat.size() > 0 )
      {
        nlon = to_int(lonlat[0]);
        nlat = to_int(lonlat[1]);
      }
    }
  }
  if ( uid == "rgg" )
  {
    if( tokens.size() > 1 ) uid = key;
    else
      uid = "rgg.N"+Translator<int,std::string>()( Resource<int>("-N",N ) );
  }
  if( !N ) N = Resource<int>("-N",N);
  if( !nlon ) nlon = Resource<int>("-nlon",nlon);
  if( !nlat ) nlat = Resource<int>("-nlat",nlat);

  if( !Factory<Grid>::instance().exists(uid) )
  {
    Log::error() << "Grid " << uid << " was not found." << std::endl;
    return ReducedGrid::Ptr();
  }

  GridSpec specs(uid);
  ReducedGrid::Ptr grid;
  if( N ) specs.set("N",N);
  if( nlon ) specs.set("nlon",nlon);
  if( nlat ) specs.set("nlat",nlat);

  grid = ReducedGrid::Ptr( ReducedGrid::create(specs) );
  return grid;
}


void AtlasGrids::run()
{
  if( !do_run ) return;

  atlas::grids::load();

  if( all )
  {
    std::vector<std::string> keys = Factory<Grid>::instance().keys();
    Log::info() << "usage: atlas-grids GRID [OPTION]... [--help]\n" << std::endl;
    Log::info() << "Available grids:" << std::endl;
    for( int i=0; i< keys.size(); ++i )
    {
      Log::info() << "  -- " << keys[i] << std::endl;
    }
  }

  if( !key.empty() )
  {
    ReducedGrid::Ptr grid = create_grid(key);

    if( !grid ) return;

    if( info )
    {
      double deg, km;
      Log::info() << "Grid " << key << std::endl;
      Log::info() << "   type:                               "
                  << grid->grid_type() << std::endl;
      Log::info() << "   uid:                                "
                  << grid->uid() << std::endl;
      if( grid->grid_type() == GaussianGrid::grid_type_str() )
      {
        Log::info() << "   N number:                           "
                    << dynamic_cast<GaussianGrid*>(grid.get())->N() << std::endl;
      }
      if( grid->grid_type() == ReducedGaussianGrid::grid_type_str() )
      {
        Log::info() << "   N number:                           "
                    << dynamic_cast<ReducedGaussianGrid*>(grid.get())->N() << std::endl;
      }
      Log::info() << "   number of points:                   "
                  << grid->npts() << std::endl;
      Log::info() << "   number of latitudes (N-S):          "
                  << grid->nlat() << std::endl;
      Log::info() << "   number of longitudes (max):         "
                  << grid->nlonmax() << std::endl;
      deg = 180./grid->nlat();
      km  = deg*40075./360.;
      Log::info() << "   approximate resolution N-S:         "
                  << std::setw(10) << std::fixed << deg << " deg   " << km << " km " << std::endl;
      deg = 360./grid->npts_per_lat()[grid->nlat()/2];
      km  = deg*40075./360.;
      Log::info() << "   approximate resolution E-W equator: "
                  << std::setw(10) << std::fixed << deg << " deg   " << km << " km " << std::endl;
      deg =  360.*std::cos(grid->lat(grid->nlat()/4)*M_PI/180.)/grid->npts_per_lat()[grid->nlat()/4];
      km  = deg*40075./360.;
      Log::info() << "   approximate resolution E-W midlat:  "
                  << std::setw(10) << std::fixed << deg << " deg   " << km << " km " << std::endl;
      deg = 360.*std::cos(grid->lat(0)*M_PI/180.)/grid->npts_per_lat()[0];
      km  = deg*40075./360.;
      Log::info() << "   approximate resolution E-W pole:    "
                  << std::setw(10) << std::fixed << deg << " deg   " << km << " km " << std::endl;
      Log::info() << "   spectral truncation -- linear:      "
                  << grid->nlat() - 1 << std::endl;
      Log::info() << "   spectral truncation -- quadratic:   "
                  << static_cast<int>(std::floor(2./3.*grid->nlat()+0.5))-1 << std::endl;
      Log::info() << "   spectral truncation -- cubic:       "
                  << static_cast<int>(std::floor(0.5*grid->nlat()+0.5))-1 << std::endl;
      Log::info() << "   bounding box -- lat N-S (deg):      "
                  << std::setw(10) << grid->bounding_box().max().lat() << ",  "
                  << std::setw(10) << grid->bounding_box().min().lat() << std::endl;
      Log::info() << "   bounding box -- lon W-E (deg):      "
                  << std::setw(10) << grid->bounding_box().min().lon() << ",  "
                  << std::setw(10) << grid->bounding_box().max().lon() << std::endl;
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
