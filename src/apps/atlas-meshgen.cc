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

#include "eckit/exception/Exceptions.h"
#include "eckit/config/Resource.h"
#include "eckit/runtime/Tool.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/parser/Tokenizer.h"
#include "eckit/geometry/Point3.h"
#include "atlas/atlas.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/mesh/actions/GenerateMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildStatistics.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/grid/grids.h"
#include "atlas/functionspace/NodeColumns.h"

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;
using namespace atlas::mesh::actions;
using namespace atlas::grid;
using namespace atlas::functionspace;
using namespace atlas::mesh;

//------------------------------------------------------------------------------------------------------

class Meshgen2Gmsh : public eckit::Tool {

  virtual void run();

public:

  Meshgen2Gmsh(int argc,char **argv): eckit::Tool(argc,argv)
  {
    bool help = Resource< bool >("--help",false);

    do_run = true;

    std::string help_str =
        "NAME\n"
        "       atlas-meshgen - Mesh generator for ReducedGrid compatible meshes\n"
        "\n"
        "SYNOPSIS\n"
        "       atlas-meshgen GRID [OPTION]... [--help] \n"
        "\n"
        "DESCRIPTION\n"
        "\n"
        "       GRID: unique identifier for grid \n"
        "           Example values: rgg.N80, rgg.TL159, gg.N40, ll.128x64\n"
        "\n"
        "       -o       Output file for mesh\n"
        "\n"
        "AUTHOR\n"
        "       Written by Willem Deconinck.\n"
        "\n"
        "ECMWF                        November 2014"
        ;
    if( help )
    {
      Log::info() << help_str << std::endl;
      do_run = false;
    }

    if( argc == 1 )
    {
      Log::info() << "usage: atlas-meshgen GRID [OPTION]... [--help]" << std::endl;
      do_run = false;
    }

    atlas_init(argc,argv);

    key = "";
    for( int i=0; i<argc; ++i )
    {
      if( i==1 && argv[i][0] != '-' )
      {
        key = std::string(argv[i]);
      }
    }

    edges      = Resource< bool> ( "--edges", false );
    stats      = Resource< bool> ( "--stats", false );
    info       = Resource< bool> ( "--info", false );
    halo       = Resource< int > ( "--halo", 0 );
    surfdim    = Resource< int > ( "--surfdim", 2 );
    brick      = Resource< int > ( "--brick", false );

    path_out = Resource<std::string> ( "-o", "" );
    if( path_out.asString().empty() && do_run )
      throw UserError(Here(),"missing output filename, parameter -o");

    if( edges )
      halo = std::max(halo,1);

  }

private:

  bool do_run;
  std::string key;
  int halo;
  bool edges;
  bool brick;
  bool stats;
  bool info;
  int surfdim;
  std::string identifier;
  std::vector<long> reg_nlon_nlat;
  std::vector<long> fgg_nlon_nlat;
  std::vector<long> rgg_nlon;
  PathName path_out;
};

//------------------------------------------------------------------------------------------------------

void Meshgen2Gmsh::run()
{
  if( !do_run ) return;
  grid::load();

  ReducedGrid::Ptr grid;
  try{ grid = ReducedGrid::Ptr( ReducedGrid::create(key) ); }
  catch( eckit::BadParameter& err ){}

  if( !grid ) return;
  mesh::Mesh::Ptr mesh;

  mesh = mesh::Mesh::Ptr( generate_mesh(*grid) );

  SharedPtr<functionspace::NodeColumns> nodes_fs( new functionspace::NodeColumns(*mesh,Halo(halo)) );
  nodes_fs->checksum(mesh->nodes().lonlat());

  Log::info() << "  checksum lonlat : " << nodes_fs->checksum(mesh->nodes().lonlat()) << std::endl;
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
  if( surfdim == 3 )
  {
    mesh::actions::BuildXYZField("xyz")(*mesh);
    gmsh.options.set("nodes",std::string("xyz"));
  }
  gmsh.write( *mesh, path_out );
  atlas_finalize();
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  Meshgen2Gmsh tool(argc,argv);
  return tool.start();
}
