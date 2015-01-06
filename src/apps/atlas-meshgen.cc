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
#include "eckit/filesystem/PathName.h"
#include "eckit/parser/Tokenizer.h"
#include "atlas/atlas.h"
#include "atlas/io/Gmsh.h"
#include "atlas/actions/GenerateMesh.h"
#include "atlas/actions/BuildEdges.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/actions/BuildHalo.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildDualMesh.h"
#include "atlas/actions/BuildStatistics.h"
#include "atlas/mpl/MPL.h"
#include "atlas/Mesh.h"
#include "atlas/grids/grids.h"
#include "atlas/GridSpec.h"

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;
using namespace atlas::actions;
using namespace atlas::grids;

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

    edges      = Resource< bool> ( "-edges", false );
    halo       = Resource< int > ( "-halo", 0 );
    surfdim    = Resource< int > ( "-surfdim", 2 );

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
  int surfdim;
  std::string identifier;
  std::vector<long> reg_nlon_nlat;
  std::vector<long> fgg_nlon_nlat;
  std::vector<long> rgg_nlon;
  PathName path_out;
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

void Meshgen2Gmsh::run()
{
  if( !do_run ) return;
  grids::load();

  ReducedGrid::Ptr grid = create_grid(key);

  if( !grid ) return;
  Mesh::Ptr mesh;

  mesh = Mesh::Ptr( generate_mesh(*grid) );

  build_nodes_parallel_fields(mesh->function_space("nodes"));
  build_periodic_boundaries(*mesh);

  if( halo )
  {
    build_halo(*mesh,halo);
    renumber_nodes_glb_idx(mesh->function_space("nodes"));
  }
  mesh->function_space("nodes").parallelise();
  ArrayView<double,2> coords( mesh->function_space("nodes").field("coordinates") );
  Log::info() << "  checksum coordinates : " << mesh->function_space("nodes").checksum().execute( coords ) << std::endl;
  if( edges )
  {
    build_edges(*mesh);
    build_pole_edges(*mesh);
    build_edges_parallel_fields(mesh->function_space("edges"),mesh->function_space("nodes"));
    build_median_dual_mesh(*mesh);
    build_statistics(*mesh);
  }

  atlas::io::Gmsh gmsh;
  gmsh.options.set("surfdim",surfdim);
  gmsh.write( *mesh, path_out );
  atlas_finalize();
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  Meshgen2Gmsh tool(argc,argv);
  tool.start();
  return 0;
}
