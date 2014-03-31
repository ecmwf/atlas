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

#include "atlas/Gmsh.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/grid/Tesselation.h"
#include "atlas/grid/GribRead.h"


//------------------------------------------------------------------------------------------------------

#if 1
#define DBG     std::cout << Here() << std::endl;
#define DBGX(x) std::cout << #x << " -> " << x << std::endl;
#else
#define DBG
#define DBGX(x)
#endif

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;

//------------------------------------------------------------------------------------------------------

class Grib2Gmsh : public eckit::Tool {

    virtual void run();

    void grib_load( const std::string& fname, atlas::Mesh& mesh, bool read_field = true );

public:

    Grib2Gmsh(int argc,char **argv): eckit::Tool(argc,argv)
    {
        in_filename = Resource<std::string>("-i","");
        if( in_filename.empty() )
            throw UserError(Here(),"missing input filename, parameter -i");

        out_filename = Resource<std::string>("-o","");
        if( out_filename.empty() )
            throw UserError(Here(),"missing output filename, parameter -o");
    }

private:

    bool gmsh;
    std::string in_filename;
    std::string out_filename;
};

//------------------------------------------------------------------------------------------------------

/// @todo this will become an expression object
void Grib2Gmsh::grib_load( const std::string& fname, atlas::Mesh& mesh, bool read_field )
{
    FILE* fh = ::fopen( fname.c_str(), "r" );
    if( fh == 0 )
        throw ReadError( std::string("error opening file ") + fname );

    int err = 0;
    grib_handle* h = grib_handle_new_from_file(0,fh,&err);

    if( h == 0 || err != 0 )
        throw ReadError( std::string("error reading grib file ") + fname );

    GribRead::read_nodes_from_grib( h, mesh );

    if( read_field )
        GribRead::read_field_from_grib( h, mesh, "field" );

    grib_handle_delete(h);

    // close file handle

    if( ::fclose(fh) == -1 )
        throw ReadError( std::string("error closing file ") + fname );
}

//------------------------------------------------------------------------------------------------------

void Grib2Gmsh::run()
{    
    std::cout.precision(std::numeric_limits< double >::digits10);
    std::cout << std::fixed;

    // input grid + field

    std::unique_ptr< atlas::Mesh > in_mesh ( new Mesh() );

    grib_load( in_filename, *in_mesh );

    atlas::Tesselation::tesselate( *in_mesh );

    atlas::Gmsh::write3dsurf( *in_mesh,out_filename );
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
    Grib2Gmsh tool(argc,argv);
    tool.start();
    return 0;
}

