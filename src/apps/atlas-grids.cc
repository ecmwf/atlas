/*
 * (C) Copyright 1996-2017 ECMWF.
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

#include "atlas/atlas.h"
#include "atlas/runtime/Log.h"
#include "atlas/grid.h"
#include "atlas/internals/AtlasTool.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/config/Resource.h"
#include "eckit/runtime/Main.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/memory/Factory.h"
#include "eckit/memory/Builder.h"
#include "eckit/log/Log.h"
#include "eckit/log/Bytes.h"
#include "eckit/parser/JSON.h"
#include "eckit/parser/Tokenizer.h"

using namespace atlas;
using namespace atlas::grid;
using eckit::JSON;
using eckit::Factory;
using eckit::SharedPtr;

//----------------------------------------------------------------------------------------------------------------------

class AtlasGrids : public AtlasTool {

  virtual bool serial() { return true; }
  virtual void execute(const Args& args);
  virtual std::string briefDescription() {
    return "Catalogue of available built-in grids";
  }
  virtual std::string usage() {
    return name() + " GRID [OPTION]... [--help,-h]";
  }
  virtual std::string longDescription() {
    return "Catalogue of available built-in grids\n"
           "\n"
           "       Browse catalogue of grids\n"
           "\n"
           "       GRID: unique identifier for grid \n"
           "           Example values: N80, F40, O24, L32\n";
  }

public:

  AtlasGrids(int argc,char **argv): AtlasTool(argc,argv)
  {
    add_option( new SimpleOption<bool>("list","List all grids. The names are possible values for the GRID argument") );

    add_option( new SimpleOption<bool>("info","List information about GRID") );

    add_option( new SimpleOption<bool>("json","Export json") );

    add_option( new SimpleOption<bool>("rtable","Export IFS rtable") );

  }

private:

  bool list;
  std::string key;
  bool info;
  bool json;
  bool rtable;
  bool do_run;
};

//------------------------------------------------------------------------------------------------------

void AtlasGrids::execute(const Args& args)
{
  key = "";
  if( args.count() ) key = args(0);

  info = false;
  args.get("info",info);

  if( info && !key.empty() ) do_run = true;

  json = false;
  args.get("json",json);
  if( json && !key.empty() ) do_run = true;

  rtable = false;
  args.get("rtable",rtable);
  if( rtable && !key.empty() ) do_run = true;

  list = false;
  args.get("list",list);
  if( list) do_run = true;

  if( !key.empty() && do_run == false )
  {
    Log::error() << "Option wrong or missing after '" << key << "'" << std::endl;
  }
  if( list )
  {
    std::vector<std::string> keys = Factory<Grid::grid_t>::instance().keys();
    Log::info() << "usage: atlas-grids GRID [OPTION]... [--help]\n" << std::endl;
    Log::info() << "Available grids:" << std::endl;
    for(size_t i = 0; i < keys.size(); ++i)
    {
      Log::info() << "  -- " << keys[i] << std::endl;
    }
  }

  if( !key.empty() )
  {
    StructuredGrid grid;
    try{ grid = Grid(key); }
    catch( eckit::BadParameter& err ){}

    if( !grid ) return;

    if( info )
    {
      double deg, km;
      Log::info() << "Grid " << key << std::endl;
//      Log::info() << "   type:                               "
//                  << g.gridType() << std::endl;
      Log::info() << "   name:                               "
                  << grid.name() << std::endl;
      Log::info() << "   uid:                                "
                  << grid.uid() << std::endl;
      // Log::info() << "   N number:                           "
      //             << grid->ny()/2 << std::endl;
      Log::info() << "   number of points:                   "
                  << grid.npts() << std::endl;
      Log::info() << "   number of latitudes (N-S):          "
                  << grid.ny() << std::endl;
      Log::info() << "   number of longitudes (max):         "
                  << grid.nxmax() << std::endl;

      deg = (grid.y().front()-grid.y().back())/(grid.ny()-1);
      km  = deg*40075./360.;
      Log::info() << "   approximate resolution N-S:         "
                  << std::setw(10) << std::fixed << deg << " deg   " << km << " km " << std::endl;


      deg = 360./static_cast<double>(grid.nx(grid.ny()/2));
      km  = deg*40075./360.;
      Log::info() << "   approximate resolution E-W equator: "
                  << std::setw(10) << std::fixed << deg << " deg   " << km << " km " << std::endl;

      deg =  360.*std::cos(grid.y(grid.ny()/4)*M_PI/180.)/
             static_cast<double>(grid.nx(grid.ny()/4));
      km  = deg*40075./360.;
      Log::info() << "   approximate resolution E-W midlat:  "
                  << std::setw(10) << std::fixed << deg << " deg   " << km << " km " << std::endl;

      deg = 360.*std::cos(grid.y().front()*M_PI/180.)/static_cast<double>(grid.nx().front());
      km  = deg*40075./360.;

      size_t memsize = grid.npts() * sizeof(double);

      Log::info() << "   memory footprint per field:                   "
                  << eckit::Bytes(memsize) << std::endl;

      Log::info() << "   approximate resolution E-W pole:    "
                  << std::setw(10) << std::fixed << deg << " deg   " << km << " km " << std::endl;

      Log::info() << "   spectral truncation -- linear:      "
                  << grid.ny() - 1 << std::endl;
      Log::info() << "   spectral truncation -- quadratic:   "
                  << static_cast<int>(std::floor(2./3.*grid.ny()+0.5))-1 << std::endl;
      Log::info() << "   spectral truncation -- cubic:       "
                  << static_cast<int>(std::floor(0.5*grid.ny()+0.5))-1 << std::endl;

    }
    if( json )
    {
      std::stringstream stream;
      JSON js(stream);
      js.precision(16);
      js << grid.spec();
      std::cout << stream.str() << std::endl;
    }

    if( rtable )
    {
      std::stringstream stream;
      stream << "&NAMRGRI\n";
      for(size_t j = 0; j < grid.ny(); ++j)
        stream << " NRGRI("<< std::setfill('0') << std::setw(5) << 1+j <<")="<< std::setfill(' ') << std::setw(5) << grid.nx(j) <<",\n";
      stream << "/" << std::flush;
      std::cout << stream.str() << std::endl;
    }
  }
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  AtlasGrids tool(argc,argv);
  return tool.start();
}
