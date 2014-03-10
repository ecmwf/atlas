#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>

#include "atlas/Gmsh.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/MeshGen.hpp"

//------------------------------------------------------------------------------------------------------

using namespace atlas;

//------------------------------------------------------------------------------------------------------

#define NLATS 256
#define NLONG 256

//------------------------------------------------------------------------------------------------------

std::vector< Point3 >* generate_ll_points( size_t nlats, size_t nlong )
{
    // generate lat/long points

    std::vector< Point3 >* pts = new std::vector< Point3 >( NLATS * NLONG );

    const double lat_inc = 180. / NLATS;
    const double lat_start = -90 + 0.5*lat_inc;
//    const double lat_end   = 90. - 0.5*lat_inc;

    const double lon_inc = 360. / NLONG;
    const double lon_start = 0.5*lon_inc;
//    const double lon_end   = 360. - 0.5*lon_inc;

    double lat = lat_start;
    double lon = lon_start;
    for( size_t ilat = 0; ilat < NLATS; ++ilat )
    {
        lon = lon_start;
        for( size_t jlon = 0; jlon < NLATS; ++jlon )
        {
//            std::cout << lat << " " << lon << std::endl;

            atlas::latlon_to_3d( lat, lon, (*pts)[ ilat*NLATS + jlon ].x );

            lon += lon_inc;
        }
        lat += lat_inc;
    }

    return pts;
}

//------------------------------------------------------------------------------------------------------

int main()
{
    std::vector< Point3 >* pts;

    pts = generate_ll_points(NLATS, NLONG);

    std::cout << "generated " << pts->size() << " points" << std::endl;

    Mesh* mesh = atlas::MeshGen::generate_from_points( *pts );

    atlas::Gmsh::write3dsurf(*mesh, std::string("earth.msh") );

    delete pts;
    delete mesh;

    return 0;
}
