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
#include "atlas/internals/AtlasTool.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/BuildStatistics.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/mesh/actions/BuildTorusXYZField.h"
#include "atlas/mesh/generators/MeshGenerator.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/util/io/Gmsh.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/geometry/Point3.h"
#include "eckit/mpi/ParallelContextBehavior.h"
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
using atlas::util::Config;

//------------------------------------------------------------------------------

class Meshgen2Gmsh : public AtlasTool {

  virtual void execute(const Args& args);
  virtual std::string briefDescription() {
    return "Mesh generator for Structured compatible meshes";
  }
  virtual std::string usage() {
    return name() + " (--grid.name=name|--grid.json=path) [OPTION]... OUTPUT [--help]";
  }

public:

  Meshgen2Gmsh(int argc,char **argv);

private:

  std::string key;
  long halo;
  bool edges;
  bool brick;
  bool stats;
  bool info;
  bool with_pole;
  bool stitch_pole;
  bool ghost;
  bool binary;
  std::string identifier;
  std::vector<long> reg_nlon_nlat;
  std::vector<long> fgg_nlon_nlat;
  std::vector<long> rgg_nlon;
  PathName path_in;
  PathName path_out;

};

//-----------------------------------------------------------------------------

Meshgen2Gmsh::Meshgen2Gmsh(int argc,char **argv): AtlasTool(argc,argv)
{
  add_option( new SimpleOption<std::string>("grid.name","Grid unique identifier\n"
    +indent()+"     Example values: N80, F40, O24, L32") );
  add_option( new SimpleOption<PathName>("grid.json","Grid described by json file") );
  add_option( new SimpleOption<double>("angle","Maximum element-edge slant deviation from meridian in degrees. \n"
    +indent()+"     Value range between 0 and 30\n"
    +indent()+"         0: Mostly triangular, with only perfect quads\n"
    +indent()+"        30: Mostly skewed quads with only triags when skewness becomes too large\n"
    +indent()+"        -1: Only triangles") );

  add_option( new SimpleOption<bool>("3d","Output mesh as sphere, and generate mesh connecting East and West in case serial") );
  add_option( new SimpleOption<bool>("include_pole","Include pole point") );
  add_option( new SimpleOption<bool>("patch_pole","Patch poles with elements.") );
  add_option( new SimpleOption<bool>("ghost","Output ghost elements") );
  add_option( new Separator("Advanced") );
  add_option( new SimpleOption<long>("halo","Halo size") );
  add_option( new SimpleOption<bool>("edges","Build edge datastructure") );
  add_option( new SimpleOption<bool>("brick","Build brick dual mesh") );
  add_option( new SimpleOption<bool>("stats","Write statistics file") );
  add_option( new SimpleOption<bool>("info","Write Info") );
  add_option( new SimpleOption<bool>("binary","Write binary file") );
  add_option( new SimpleOption<std::string>("generator","Mesh generator") );
  add_option( new SimpleOption<std::string>("partitioner","Mesh partitioner") );
  add_option( new SimpleOption<bool>("periodic_x","periodic mesh in x-direction") );
  add_option( new SimpleOption<bool>("periodic_y","periodic mesh in y-direction") );
  add_option( new SimpleOption<bool>("torus","Output mesh as torus") );
}

//-----------------------------------------------------------------------------

void Meshgen2Gmsh::execute(const Args& args)
{
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
  brick = false;
  args.get("brick",brick);
  ghost = false;
  args.get("ghost",ghost);
  binary = false;
  args.get("binary",binary);

  std::string path_in_str = "";
  if( args.get("grid.json",path_in_str) ) path_in = path_in_str;

  if( args.count() )
    path_out = args(0);
  else
    path_out = "mesh.msh";

  if( path_in_str.empty() && key.empty() ) {
    Log::warning() << "missing argument --grid.name or --grid.json" << std::endl;
    Log::warning() << "Usage: " << usage() << std::endl;
    return;
  }


  if( edges )
    halo = std::max(halo,1l);

  std::string meshgenerator_type("Structured");
  args.get("generator",meshgenerator_type);
  eckit::LocalConfiguration meshgenerator_config( args );
  if( eckit::mpi::size() > 1 || edges )
    meshgenerator_config.set("3d",false);

  SharedPtr<Structured> grid;
  if( key.size() )
  {
    try{ grid.reset( Structured::create(key) ); }
    catch( eckit::BadParameter& e ){}
  }
  else if( path_in.path().size() )
  {
    Log::info() << "Creating grid from file " << path_in << std::endl;
    Log::debug() << atlas::util::Config(path_in) << std::endl;
    try{ grid.reset( Structured::create( atlas::util::Config(path_in) ) ); }
    catch( eckit::BadParameter& e ){}
  }
  else
  {
    Log::error() << "No grid specified." << std::endl;
  }

  if( !grid ) return;
  
  Log::debug() << "Domain spec: " << grid->domain().spec() << std::endl;
  Log::debug() << "Domain global: " << std::string(grid->domain().isGlobal() ? "yes" : "no") << std::endl;
  
  SharedPtr<mesh::generators::MeshGenerator> meshgenerator (
      mesh::generators::MeshGenerator::create(meshgenerator_type,meshgenerator_config) );


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

  if( edges )
  {
    build_edges(*mesh);
    build_pole_edges(*mesh);
    build_edges_parallel_fields(*mesh);
    if( brick )
      build_brick_dual_mesh(*grid, *mesh);
    else
      build_median_dual_mesh(*mesh);
  }

  if( stats )
    build_statistics(*mesh);

  bool torus=false;
  args.get("torus",torus);
  if( torus ) {
    dim_3d = true;
    Log::debug() << "Building xyz representation for nodes on torus" << std::endl;
    mesh::actions::BuildTorusXYZField("xyz")(*mesh,&grid->domain(),5.,2.,grid->nlonmax(),grid->nlat());
  }

  atlas::output::Gmsh gmsh( path_out , Config
    ("info",info)
    ("ghost",ghost)
    ("coordinates", dim_3d ? "xyz" : "geolonlat" )
    ("edges",edges )
    ("binary",binary )
  );
  Log::info() << "Writing mesh to gmsh file \"" << path_out << "\" generated from grid \"" << grid->shortName() << "\"" << std::endl;
  gmsh.write( *mesh );
}

//------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  Meshgen2Gmsh tool(argc,argv);
  return tool.start();
}
