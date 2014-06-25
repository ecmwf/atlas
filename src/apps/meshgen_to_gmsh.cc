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
#include "eckit/value/ListContent.h"

#include "atlas/io/Gmsh.hpp"
#include "atlas/actions/GenerateMesh.hpp"
#include "atlas/mpl/MPL.hpp"
#include "atlas/mesh/Mesh.hpp"

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
      lon  = Resource<int>("-lon",-1);
      lat  = Resource<int>("-lat",-1);

      identifier = Resource<std::string>("-rgg","");
      if( identifier.empty() && (lon<0 || lat<0) )
          throw UserError(Here(),"missing input mesh identifier, parameter -rgg");

      path_out = Resource<std::string>("-o","");
      if( path_out.asString().empty() )
        throw UserError(Here(),"missing output filename, parameter -o");

    }

private:

    std::string identifier;
    int lon, lat;
    PathName path_out;
};

//------------------------------------------------------------------------------------------------------

void Meshgen2Gmsh::run()
{
    MPL::init();
    Mesh::Ptr mesh;
    if( !identifier.empty() )
      mesh = Mesh::Ptr( generate_reduced_gaussian_grid(identifier) );
    else
      mesh = Mesh::Ptr( generate_latlon_grid(lon,lat) );
    atlas::Gmsh::write( *mesh, path_out );
    MPL::finalize();
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
    Meshgen2Gmsh tool(argc,argv);
    tool.start();
    return 0;
}

