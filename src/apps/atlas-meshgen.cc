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
#include "atlas/atlas.h"
#include "atlas/io/Gmsh.h"
#include "atlas/actions/GenerateMesh.h"
#include "atlas/actions/BuildEdges.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/actions/BuildHalo.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildDualMesh.h"
#include "atlas/mpl/MPL.h"
#include "atlas/Mesh.h"
#include "atlas/grids/grids.h"

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;
using namespace atlas::actions;

//------------------------------------------------------------------------------------------------------

class Meshgen2Gmsh : public eckit::Tool {

  virtual void run();

public:

  Meshgen2Gmsh(int argc,char **argv): eckit::Tool(argc,argv)
  {
    atlas_init(argc,argv);

    rgg_nlon   = Resource< std::vector<long> > ( "-rgg_nlon", std::vector<long>() );
    reg_nlon_nlat  = Resource< std::vector<long> > ( "-reg", std::vector<long>() );
    fgg_nlon_nlat  = Resource< std::vector<long> > ( "-fgg", std::vector<long>() );
    identifier = Resource< std::string       > ( "-rgg", "" );
    edges      = Resource< bool > ( "-edges", false );
    halo       = Resource< int > ( "-halo", 0 );
    surfdim    = Resource< int > ( "-surfdim", 2 );

    int col = Resource< int >( "-col", 0);
    if( col )
    {
      rgg_nlon.resize(col);
      for( int jlat=0; jlat<(col-4); ++jlat )
      {
        rgg_nlon[jlat] = 20+4*jlat;
      }
      for( int jlat=(col-4); jlat<col; ++jlat )
      {
        rgg_nlon[jlat] = rgg_nlon[(col-4)-1];
      }
    }

    if( identifier.empty() && reg_nlon_nlat.empty() && fgg_nlon_nlat.empty() && rgg_nlon.empty() )
    {
      throw UserError(Here(),"missing input mesh identifier, parameter -rgg or -reg");
    }
    path_out = Resource<std::string> ( "-o", "" );
    if( path_out.asString().empty() )
      throw UserError(Here(),"missing output filename, parameter -o");

    if( edges )
      halo = std::max(halo,1);

  }

private:

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

void Meshgen2Gmsh::run()
{
  Mesh::Ptr mesh;
  if( !identifier.empty() )
    mesh = Mesh::Ptr( generate_reduced_gaussian_grid(identifier) );
  else if( !fgg_nlon_nlat.empty() )
    mesh = Mesh::Ptr( generate_full_gaussian_grid(fgg_nlon_nlat[0],fgg_nlon_nlat[1]) );
  else if( !reg_nlon_nlat.empty() )
    mesh = Mesh::Ptr( generate_regular_grid(reg_nlon_nlat[0],reg_nlon_nlat[1]) );
  else if ( !rgg_nlon.empty() )
    mesh = Mesh::Ptr( generate_reduced_gaussian_grid(rgg_nlon) );
  else
    throw UserError(Here(),"Could not generate mesh with given input");

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
