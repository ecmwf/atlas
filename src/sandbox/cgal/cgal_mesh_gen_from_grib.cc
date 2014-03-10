#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>

#include "grib_api.h"

#include "atlas/Gmsh.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/MeshGen.hpp"

//------------------------------------------------------------------------------------------------------

using namespace atlas;

//------------------------------------------------------------------------------------------------------

#define NLATS 256
#define NLONG 256

//------------------------------------------------------------------------------------------------------

std::vector< Point3 >* read_ll_points_from_grib( const std::string& filename )
{
    int err = 0;

    // points to read

    std::vector< Point3 >* pts = new std::vector< Point3 >();

    // load grib file

    FILE* f = ::fopen( filename.c_str(), "r" );
    if( f == 0 )
        throw std::string("error opening file");

    grib_handle* h = grib_handle_new_from_file(0,f,&err);

    if( h == 0 || err != 0 )
        throw std::string("error reading grib");

    grib_iterator *i = grib_iterator_new(h, 0, &err);

    double lat   = 0.;
    double lon   = 0.;
    double value = 0.;

    while(grib_iterator_next(i,&lat,&lon,&value))
    {
        while(lon < 0)    lon += 360;
        while(lon >= 360) lon -= 360;

        pts->push_back( Point3() );

        atlas::latlon_to_3d( lat, lon, pts->back().x );
    }
    grib_iterator_delete(i);

    if( ::fclose(f) == -1 )
        throw std::string("error closing file");

    return pts;
}

//------------------------------------------------------------------------------------------------------

int main()
{
    std::vector< Point3 >* pts;

    pts = read_ll_points_from_grib( "data.grib" );

    std::cout << "generated " << pts->size() << " points" << std::endl;

    Mesh* mesh = atlas::MeshGen::generate_from_points( *pts );

    atlas::Gmsh::write3dsurf(*mesh, std::string("earth.msh") );

    delete pts;
    delete mesh;

    return 0;
}
