#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>

#include "atlas/Parameters.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/Field.hpp"
#include "atlas/FunctionSpace.hpp"

#include "atlas/Gmsh.hpp"

//------------------------------------------------------------------------------------------------------

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Segment_3                              Segment_3;

// define point creator

typedef K::Point_3                                Point_3;
typedef CGAL::Creator_uniform_3<double, Point_3>  PointCreator;

//------------------------------------------------------------------------------------------------------

using namespace atlas;

//------------------------------------------------------------------------------------------------------

const double earthRadius = 6367.47; // from ECMWF model ...

class LL3D {
public:

    double x_[3];

public:

    LL3D(){}

    LL3D( double lat, double lon )
    {
        assign(lat,lon);
    }

    void assign(double lat, double lon)
    {
        // See http://en.wikipedia.org/wiki/Geodetic_system#From_geodetic_to_ECEF
        double& X = x_[0];
        double& Y = x_[1];
        double& Z = x_[2];

        double h = 0; // Altitude

        double a = earthRadius; // 6378137.0 ;       //  WGS84 semi-major axis
        double e2 = 0;          // 6.69437999014E-3; // WGS84 first numerical eccentricity sqared

        double phi = lat / 180.0 * M_PI;
        double lambda = lon / 180.0 * M_PI;

        double cos_phi = cos(phi);
        double sin_phi = sin(phi);
        double cos_lambda = cos(lambda);
        double sin_lambda = sin(lambda);

        double N_phi = a/sqrt(1-e2*sin_phi*sin_phi);

        X = (N_phi + h) * cos_phi * cos_lambda;
        Y = (N_phi + h) * cos_phi * sin_lambda;
        Z = (N_phi * (1-e2) + h) * sin_phi;
    }

//    friend std::ostream& operator<<(std::ostream& s,const LLPoint2& p)
//    {
//        s << '(' << p.lat_ << "," << p.lon_ << ')';
//        return s;
//    }

};

#define NLATS 25
#define NLONG 25

//------------------------------------------------------------------------------------------------------

std::vector< LL3D >* generate_ll_points( size_t nlats, size_t nlong )
{
    // generate lat/long points

    std::vector< LL3D >* pts = new std::vector< LL3D >( NLATS * NLONG );

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

            (*pts)[ ilat*NLATS + jlon ].assign( lat, lon );
            lon += lon_inc;
        }
        lat += lat_inc;
    }

    return pts;
}

//------------------------------------------------------------------------------------------------------

Polyhedron_3* create_convex_hull_from_points( const std::vector< LL3D >& pts )
{
    Polyhedron_3* poly = new Polyhedron_3();

    // insertion from a vector :

    std::vector<Point_3> vertices( pts.size() );
    for( size_t i = 0; i < vertices.size(); ++i )
    {
        vertices[i] = Point_3( pts[i].x_[0], pts[i].x_[1], pts[i].x_[2] );
    }

    // compute convex hull of non-collinear points

    CGAL::convex_hull_3( vertices.begin(), vertices.end(), *poly );

    return poly;
}

//------------------------------------------------------------------------------------------------------

int main()
{
    std::vector< LL3D >* pts;

    pts = ( generate_ll_points(NLATS, NLONG) );

    std::cout << "generated " << pts->size() << " points" << std::endl;

    // define polyhedron to hold convex hull

    Polyhedron_3* poly = create_convex_hull_from_points( *pts );

    std::cout << "convex hull " << poly->size_of_vertices() << " vertices" << std::endl;

    assert( poly->size_of_vertices() == pts->size() );

    delete pts;
    delete poly;

    return 0;
}
