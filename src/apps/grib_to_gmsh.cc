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

#include "atlas/util/io/Gmsh.h"
#include "atlas/mesh/Mesh.h"
//#include "atlas/grid/Tesselation.h"
//#include "atlas/grid/GribRead.h"
#include "atlas/field/FieldSet.h"

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;
using namespace atlas::grid;

//------------------------------------------------------------------------------------------------------

class Grib2Gmsh : public eckit::Tool {

    virtual void run();

public:

    Grib2Gmsh(int argc,char **argv): eckit::Tool(argc,argv)
    {
        path_in = Resource<std::string>("-i","");
        if( path_in.asString().empty() )
            throw UserError(Here(),"missing input filename, parameter -i");

        path_out = Resource<std::string>("-o","");
        if( path_out.asString().empty() )
            throw UserError(Here(),"missing output filename, parameter -o");
    }

private:

    PathName path_in;
    PathName path_out;
};

//------------------------------------------------------------------------------------------------------

void Grib2Gmsh::run()
{    
    std::cout.precision(std::numeric_limits< double >::digits10);
    std::cout << std::fixed;

    field::FieldSet::Ptr fs_inp;

    fs_inp.reset( new field::FieldSet( path_in ) );
    if( fs_inp->empty() )
        throw UserError("Input fieldset is empty", Here());

    Tesselation::tesselate( fs_inp->grid() );

    atlas::util::io::Gmsh::write3dsurf( fs_inp->grid().mesh(), path_out );
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
    Grib2Gmsh tool(argc,argv);
    tool.start();
    return 0;
}

