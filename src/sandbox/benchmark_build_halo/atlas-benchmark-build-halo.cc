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
#include "eckit/log/Timer.h"

#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/mesh.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/output/detail/GmshIO.h"
#include "atlas/mesh/actions/BuildHalo.h"

//------------------------------------------------------------------------------

using namespace atlas;
using namespace atlas::grid;
using atlas::util::Config;
using eckit::PathName;

struct TimerStats
{
  TimerStats(const std::string& _name = "timer")
  {
    max = -1;
    min = -1;
    avg = 0;
    cnt = 0;
    name = _name;
  }
  void update(eckit::Timer& timer)
  {
    double t = timer.elapsed();
    if( min < 0 ) min = t;
    if( max < 0 ) max = t;
    min = std::min(min, t);
    max = std::max(max, t);
    avg = (avg*cnt+t)/(cnt+1);
    ++cnt;
  }
  std::string str()
  {
    std::stringstream stream;
    stream << name << ": min, max, avg -- " << min << ", " << max << ", " << avg;
    return stream.str();
  }
  std::string name;
  double max;
  double min;
  double avg;
  int cnt;
};

//------------------------------------------------------------------------------

class Tool : public AtlasTool {

  virtual void execute(const Args& args);
  virtual std::string briefDescription() {
    return "Tool to generate a python script that plots the grid-distribution of a given grid";
  }
  virtual std::string usage() {
    return name() + " (--grid.name=name|--grid.json=path) [OPTION]... OUTPUT [--help]";
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
  add_option( new SimpleOption<std::string>("grid.name","Grid unique identifier\n"
    +indent()+"     Example values: N80, F40, O24, L32") );
  add_option( new SimpleOption<PathName>("grid.json","Grid described by json file") );
  add_option( new SimpleOption<long>("halo","Number of halos") );
}

//-----------------------------------------------------------------------------

void Tool::execute(const Args& args)
{
  key = "";
  args.get("grid.name",key);

  std::string path_in_str = "";
  if( args.get("grid.json",path_in_str) ) path_in = path_in_str;

  StructuredGrid grid;
  if( key.size() )
  {
    try{ grid = Grid(key); }
    catch( eckit::BadParameter& e ){}
  }
  else if( path_in.path().size() )
  {
    Log::info() << "Creating grid from file " << path_in << std::endl;
    Log::debug() << Config(path_in) << std::endl;
    try{ grid = Grid( Config(path_in) ); }
    catch( eckit::BadParameter& e ){}
  }
  else
  {
    Log::error() << "No grid specified." << std::endl;
  }

  if( !grid ) return;

  Log::debug() << "Domain: "   << grid.domain() << std::endl;
  Log::debug() << "Periodic: " << grid.periodic() << std::endl;

  MeshGenerator meshgenerator("structured");

  size_t halo = args.getLong("halo",1);

  size_t iterations = 10;
  parallel::mpi::comm().barrier();
  TimerStats timer_stats;
  for( size_t i=0; i<iterations; ++i )
  {
    Mesh mesh = meshgenerator.generate(grid);
    parallel::mpi::comm().barrier();
    eckit::Timer timer;
    mesh::actions::build_halo( mesh, halo );
    parallel::mpi::comm().barrier();
    timer.stop();
    timer_stats.update(timer);
    Log::info() << "iteration " << std::setw(2) << i << " : " << std::setprecision(5) << std::fixed << timer.elapsed() << " seconds"<< std::endl;
  }
  Log::info() << "Timer Statistics:\n"
              << "  min: " << std::setprecision(5) << std::fixed << timer_stats.min
              << "  max: " << std::setprecision(5) << std::fixed << timer_stats.max
              << "  avg: " << std::setprecision(5) << std::fixed << timer_stats.avg << std::endl;

}

//------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  Tool tool(argc,argv);
  return tool.start();
}
