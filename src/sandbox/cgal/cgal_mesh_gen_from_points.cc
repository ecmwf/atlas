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

typedef K::Vector_3                               Vector_3;
typedef K::FT                                     FT;
typedef K::Segment_3                              Segment_3;
typedef K::Point_3                                Point_3;

const Point_3 origin = Point_3(CGAL::ORIGIN);

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

#define NLATS 100
#define NLONG 100

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

atlas::Mesh* cgal_polyhedron_to_atlas_mesh( Polyhedron_3& poly )
{
    bool ensure_outward_normals = true;

    Mesh* mesh = new Mesh();

    /* nodes */

    const size_t nb_nodes = poly.size_of_vertices();

    std::cout << "nb_nodes = " << nb_nodes << std::endl;

    std::vector<int> bounds(2);
    bounds[0] = Field::UNDEF_VARS;
    bounds[1] = nb_nodes;

    FunctionSpace& nodes = mesh->add_function_space( new FunctionSpace( "nodes", "Lagrange_P0", bounds ) );

    nodes.metadata().set("type",static_cast<int>(Entity::NODES));

    FieldT<double>& coords = nodes.create_field<double>("coordinates",3);

    std::map< Polyhedron_3::Vertex_const_handle, size_t > vidx;

    size_t inode = 0;
    for( Polyhedron_3::Vertex_const_iterator v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
    {
        vidx[v] = inode;

        const Polyhedron_3::Point_3& p = v->point();

        coords(XX,inode) = p.x();
        coords(YY,inode) = p.y();
        coords(ZZ,inode) = p.z();

        ++inode;
    }

    assert( inode == nb_nodes );

    /* triangles */

    const size_t nb_triags = poly.size_of_facets();

    std::cout << "nb_triags = " << nb_triags << std::endl;

    bounds[1] = nb_triags;

    FunctionSpace& triags  = mesh->add_function_space( new FunctionSpace( "triags", "Lagrange_P1", bounds ) );
    triags.metadata().set("type",static_cast<int>(Entity::ELEMS));

    FieldT<int>& triag_nodes = triags.create_field<int>("nodes",3);

    int idx[3];
    Polyhedron_3::Vertex_const_handle vts[3];

    size_t tidx = 0;
    for( Polyhedron_3::Facet_const_iterator f = poly.facets_begin(); f != poly.facets_end(); ++f )
    {

        // loop  over half-edges and take each vertex()

        size_t iedge = 0;
        Polyhedron_3::Halfedge_around_facet_const_circulator edge = f->facet_begin();

        do
        {
            Polyhedron_3::Vertex_const_handle vh = edge->vertex();
            const Polyhedron_3::Point_3& p = vh->point();

            idx[iedge] = vidx[vh];
            vts[iedge] = vh;

            ++iedge;
            ++edge;
        }
        while ( edge != f->facet_begin() && iedge < 3 );

        assert( iedge == 3 );

        if( ensure_outward_normals ) /* ensure outward pointing normal */
        {
            Vector_3 p0 ( origin, vts[0]->point() );
            Vector_3 n  = CGAL::normal( vts[0]->point(), vts[1]->point(), vts[2]->point() );

            FT innerp = n * p0;

            if( innerp < 0 ) // need to swap an edge of the triag
                std::swap( vts[1], vts[2] );
        }

        /* define the triag */

        triag_nodes(0,tidx) = F_IDX(idx[0]);
        triag_nodes(1,tidx) = F_IDX(idx[1]);
        triag_nodes(2,tidx) = F_IDX(idx[2]);

        ++tidx;

    }

    assert( tidx == nb_triags );

    return mesh;
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

    Mesh* mesh = cgal_polyhedron_to_atlas_mesh( *poly );

    atlas::Gmsh::write3dsurf(*mesh, std::string("earth.msh") );

    delete pts;
    delete poly;
    delete mesh;

    return 0;
}
