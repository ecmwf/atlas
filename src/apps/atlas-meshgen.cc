/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/grids.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/BuildStatistics.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/mesh/generators/MeshGenerator.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/util/io/Gmsh.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/geometry/Point3.h"
#include "eckit/mpi/ParallelContextBehavior.h"
#include "eckit/option/CmdArgs.h"
#include "eckit/option/Separator.h"
#include "eckit/option/SimpleOption.h"
#include "eckit/parser/Tokenizer.h"
#include "eckit/runtime/Context.h"
#include "eckit/runtime/Tool.h"


//------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;
using namespace atlas::mesh::actions;
using namespace atlas::grid;
using namespace atlas::functionspace;
using namespace atlas::mesh;
using namespace eckit::option;

//------------------------------------------------------------------------------

void usage(const std::string& tool)
{
  Log::info() << "usage: " << tool << " (--grid.name=name|--grid.json=path) [OPTION]... OUTPUT [--help]" << std::endl;
}

class Meshgen2Gmsh : public eckit::Tool {
  typedef std::vector<eckit::option::Option *> Options;

  virtual void run();

public:

  Meshgen2Gmsh(int argc,char **argv): eckit::Tool(argc,argv)
  {
    Context::instance().behavior( new eckit::mpi::ParallelContextBehavior() );
    
    if( argc == 1 )
    {
      usage( name() );
      do_run = false;
      return;
    }

    Options options;
    options.push_back( new SimpleOption<bool>("help","Print this help") );
    options.push_back( new SimpleOption<std::string>("grid.name","Grid unique identifier\n"
      "           Example values: N80, F40, O24, L32") );
    options.push_back( new SimpleOption<PathName>("grid.json","Grid described by json file") );
    options.push_back( new SimpleOption<double>("angle","Maximum element-edge slant deviation from meridian in degrees. \n"
      "           Value range between 0 and 30\n"
      "               0: Mostly triangular, with only perfect quads\n"
      "              30: Mostly skewed quads with only triags when skewness becomes too large\n"
      "              -1: Only triangles") );

    options.push_back( new SimpleOption<bool>("3d","Output mesh as sphere, and generate mesh connecting East and West in case serial") );
    options.push_back( new SimpleOption<bool>("include_pole","Include pole point") );
    options.push_back( new SimpleOption<bool>("patch_pole","Patch poles with elements.") );
    options.push_back( new SimpleOption<bool>("ghost","Output ghost elements") );
    options.push_back( new Separator("Advanced") );
    options.push_back( new SimpleOption<long>("halo","Halo size") );
    options.push_back( new SimpleOption<bool>("edges","Build edge datastructure") );
    options.push_back( new SimpleOption<bool>("brick","Build brick dual mesh") );
    options.push_back( new SimpleOption<bool>("stats","Write statistics file") );
    options.push_back( new SimpleOption<bool>("info","Write Info") );
    
    std::stringstream options_str;
    for( Options::const_iterator it = options.begin(); it!= options.end(); ++it ) {
      options_str << "    " << (**it) << "\n\n";
    }
    
    CmdArgs args(&usage, options);

    do_run = true;

    std::string help_str =
        "NAME\n"
        "       atlas-meshgen - Mesh generator for Structured compatible meshes\n"
        "\n"
        "SYNOPSIS\n"
        "       atlas-meshgen (--grid.name=name|--grid.json=path) [OPTION]... OUTPUT [--help] \n"
        "\n"
        "DESCRIPTION\n"
        "\n"
          +options_str.str()+
        "\n"
        "AUTHOR\n"
        "       Written by Willem Deconinck.\n"
        "\n"
        "ECMWF                        November 2014"
        ;
    
    bool help=false;
    args.get("help",help);
    if( help )
    {
      if( eckit::mpi::rank() == 0 )
        Log::info() << help_str << std::endl;
      do_run = false;
      return;
    }

    atlas_init(argc,argv);


    key = "";
    args.get("grid.name",key);

    edges = false;
    args.get("edges",edges);
    stats = false;
    args.get("stats",stats);
    info = false;
    args.get("info",info);
    halo       = 0;
    args.get("halo",halo);
    bool dim_3d=false;
    args.get("3d",dim_3d);
    surfdim    = dim_3d ? 3 : 2;
    brick = false;
    args.get("brick",brick);
    ghost = false;
    args.get("ghost",ghost);

    std::string path_in_str = "";
    if( args.get("grid.json",path_in_str) ) path_in = path_in_str;
      
    if( args.count() )
      path_out = args(0);
    else
      path_out = "mesh.msh";

    if( path_in.asString().empty() && key.empty() && do_run )
      Log::info() << "missing argument --grid.name or --grid.json" << std::endl;

    if( edges )
      halo = std::max(halo,1l);

    meshgenerator_config = args.get();
    if( eckit::mpi::size() > 1 )
      meshgenerator_config.set("3d",false);
  }

private:

  bool do_run;
  std::string key;
  long halo;
  bool edges;
  bool brick;
  bool stats;
  bool info;
  int surfdim;
  bool with_pole;
  bool stitch_pole;
  bool ghost;
  std::string identifier;
  std::vector<long> reg_nlon_nlat;
  std::vector<long> fgg_nlon_nlat;
  std::vector<long> rgg_nlon;
  PathName path_in;
  PathName path_out;
  
  eckit::LocalConfiguration meshgenerator_config;
  
};

//------------------------------------------------------------------------------

void Meshgen2Gmsh::run()
{
  if( !do_run ) return;
  grid::load();

  SharedPtr<global::Structured> grid;
  if( key.size() )
  {
    try{ grid.reset( global::Structured::create(key) ); }
    catch( eckit::BadParameter& e ){}
  }
  else if( path_in.path().size() )
  {
    Log::info() << "Creating grid from file " << path_in << std::endl;
    try{ grid.reset( global::Structured::create( atlas::util::Config(path_in) ) ); }
    catch( eckit::BadParameter& e ){}
  }
  else
  {
    Log::error() << "No grid specified." << std::endl;
  }

  if( !grid ) return;
  SharedPtr<mesh::generators::MeshGenerator> meshgenerator (
      mesh::generators::MeshGenerator::create("Structured",meshgenerator_config) );
  
  
  SharedPtr<mesh::Mesh> mesh;
  try {
  mesh.reset( meshgenerator->generate(*grid) );
  }
  catch ( eckit::BadParameter& e)
  {
    Log::error() << e.what() << std::endl;
    Log::error() << e.callStack() << std::endl;
    throw e;
  }
  SharedPtr<functionspace::NodeColumns> nodes_fs( new functionspace::NodeColumns(*mesh,Halo(halo)) );
  // nodes_fs->checksum(mesh->nodes().lonlat());
  // Log::info() << "  checksum lonlat : " << nodes_fs->checksum(mesh->nodes().lonlat()) << std::endl;
  
  if( edges )
  {
    build_edges(*mesh);
    build_pole_edges(*mesh);
    build_edges_parallel_fields(*mesh);
    if( brick )
      build_brick_dual_mesh(*mesh);
    else
      build_median_dual_mesh(*mesh);
  }

  if( stats )
    build_statistics(*mesh);

  atlas::util::io::Gmsh gmsh;
  gmsh.options.set("info",info);
  gmsh.options.set("ghost",ghost);
  if( surfdim == 3 )
  {
    mesh::actions::BuildXYZField("xyz")(*mesh);
    gmsh.options.set("nodes",std::string("xyz"));
  }
  Log::info() << "Writing mesh to gmsh file \"" << path_out << "\" generated from grid \"" << grid->shortName() << "\"" << std::endl;
  gmsh.write( *mesh, path_out );
  atlas_finalize();
}

//------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  Meshgen2Gmsh tool(argc,argv);
  return tool.start();
}
