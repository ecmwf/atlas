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
#include "eckit/memory/ScopedPtr.h"
#include "eckit/config/Resource.h"
#include "eckit/runtime/Tool.h"

#include "atlas/atlas.h"
#include "atlas/io/Gmsh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/Nodes.h"
#include "atlas/mesh/Field.h"
#include "atlas/mesh/FunctionSpace.h"

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;

//------------------------------------------------------------------------------------------------------

class GmshLs : public eckit::Tool {

    virtual void run();

public:

    GmshLs(int argc,char **argv): eckit::Tool(argc,argv)
    {
        atlas_init(argc,argv);
        in_filename = Resource<std::string>("-i","");
        if( in_filename.empty() )
            throw UserError(Here(),"missing input filename, parameter -i");
    }

private:

    std::string in_filename;
};

//------------------------------------------------------------------------------------------------------

void GmshLs::run()
{
    std::cout.precision(std::numeric_limits< double >::digits10);
    std::cout << std::fixed;

    // input grid + field

    atlas::Mesh::Ptr in_mesh ( atlas::Gmsh::read( in_filename ) );

    Mesh& mesh = *in_mesh;

    Nodes& nodes   = mesh.nodes();
    FieldT<double>& lonlat = nodes.field<double>( "coordinates" );
    FieldT<int>& glb_idx   = nodes.field<int>( "glb_idx" );

    size_t nb_nodes = nodes.shape(0);

    Log::info() << "nb_nodes = " << nb_nodes << std::endl;

    if( mesh.has_function_space("triags") ) Log::info() << "nb_triags = " << mesh.function_space( "triags" ).shape(0) << std::endl;
    if( mesh.has_function_space("quads") )  Log::info() << "nb_quads = "  << mesh.function_space( "quads" ).shape(0) << std::endl;
    if( mesh.has_function_space("edges") )  Log::info() << "nb_edges = "  << mesh.function_space( "edges" ).shape(0) << std::endl;

    atlas_finalize();

}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
    GmshLs tool(argc,argv);
    tool.start();
    return 0;
}

