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
#include <iomanip>

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"

#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/mesh.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/output/detail/GmshIO.h"
#include "atlas/mesh/actions/BuildHalo.h"

//------------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::grid;
using atlas::util::Config;
using eckit::PathName;

//------------------------------------------------------------------------------

class Tool : public AtlasTool {

  virtual void execute(const Args& args);
  virtual std::string briefDescription() {
    return "Tool to generate a python script that plots the grid-distribution of a given grid";
  }
  virtual std::string usage() {
    return name() + " --grid=name [OPTION]... OUTPUT [--help]";
  }

public:

  Tool(int argc,char **argv);

private:

  std::string key;
  PathName path_in;
  PathName path_out;

};

//-----------------------------------------------------------------------------

Tool::Tool(int argc,char **argv): AtlasTool(argc,argv)
{
  add_option( new SimpleOption<std::string>("grid","Grid unique identifier\n"
    +indent()+"     Example values: N80, F40, O24, L32") );
  add_option( new SimpleOption<long>("halo","Number of halos (default=1") );
}

//-----------------------------------------------------------------------------

void Tool::execute(const Args& args)
{
  Trace timer(Here(), displayName() );
  key = "";
  args.get("grid",key);

  std::string path_in_str = "";
  if( args.get("grid",path_in_str) ) path_in = path_in_str;

  StructuredGrid grid;
  if( key.size() )
  {
    try{ grid = Grid(key); }
    catch( eckit::BadParameter& e ){}
  }
  else
  {
    Log::error() << "No grid specified." << std::endl;
  }


  if( !grid ) return;

  Log::debug() << "Domain: "   << grid.domain() << std::endl;
  Log::debug() << "Periodic: " << grid.periodic() << std::endl;

  MeshGenerator meshgenerator("structured", util::Config("partitioner","equal_regions"));

  size_t halo = args.getLong("halo",1);

  size_t iterations = 10;
  for( size_t i=0; i<iterations; ++i )
  {
    ATLAS_TRACE("iteration");
    Mesh mesh = meshgenerator.generate(grid);
    mesh::actions::build_halo( mesh, halo );
    parallel::mpi::comm().barrier();
  }
  timer.stop();
  Log::info() << Trace::report() << std::endl;

}

//------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  Tool tool(argc,argv);
  return tool.start();
}
