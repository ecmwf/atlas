/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <iostream>

#include "eckit/exception/Exceptions.h"
#include "eckit/config/Resource.h"
#include "eckit/runtime/Tool.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/atlas.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/grid/grids.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/actions/GenerateMesh.h"
#include "atlas/mesh/actions/WriteLoadBalanceReport.h"

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;
using namespace atlas::mesh::actions;
using namespace atlas::grid;
using namespace atlas::functionspace;
using namespace atlas::mesh;

//------------------------------------------------------------------------------------------------------

class AtlasLoadbalance : public eckit::Tool {

  virtual void run();

public:

  AtlasLoadbalance(int argc,char **argv): eckit::Tool(argc,argv)
  {
    bool help = Resource< bool >("--help",false);

    do_run = true;

    std::string help_str =
        "NAME\n"
        "       atlas-loadbalance - <TODO>\n"
        "\n"
        "SYNOPSIS\n"
        "       atlas-loadbalance GRID [OPTION]... [--help] \n"
        "\n"
        "DESCRIPTION\n"
        "\n"
        "       GRID: unique identifier for grid \n"
        "           Example values: rgg.N80, rgg.TL159, gg.N40, ll.128x64\n"
        "\n"
        "       --halo       Output file for mesh\n"
        "\n"
        "AUTHOR\n"
        "       Written by Willem Deconinck.\n"
        "\n"
        "ECMWF                        September 2015"
        ;
    if( help )
    {
      Log::info() << help_str << std::endl;
      do_run = false;
    }

    if( argc == 1 )
    {
      Log::info() << "usage: atlas-loadbalance GRID [OPTION]... [--help]" << std::endl;
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

    halo       = Resource< int > ( "--halo", 1 );
    output     = Resource<std::string> ("--output","");
  }

private:

  bool do_run;
  std::string key;
  int halo;
  std::string output;
  std::string identifier;
};

//------------------------------------------------------------------------------------------------------

void AtlasLoadbalance::run()
{
  if( !do_run ) return;
  grid::load();

  ReducedGrid::Ptr grid;
  try{ grid = ReducedGrid::Ptr( ReducedGrid::create(key) ); }
  catch( eckit::BadParameter& err ){}

  if( !grid ) return;
  SharedPtr<mesh::Mesh> mesh( generate_mesh(*grid) );
  SharedPtr<functionspace::NodeColumns> nodes( new functionspace::NodeColumns(*mesh,Halo(halo)) );


  if( output.size() )
  {
    write_load_balance_report(*mesh,output);
  }
  else
  {
    std::stringstream s;
    write_load_balance_report(*mesh,s);

    if( eckit::mpi::rank() == 0 )
    {
      std::cout << s.str() << std::endl;
    }
  }
  atlas_finalize();
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  AtlasLoadbalance tool(argc,argv);
  return tool.start();
}
