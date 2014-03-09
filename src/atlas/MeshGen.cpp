// (C) Copyright 1996-2014 ECMWF.

#include <cmath>

#include "atlas/atlas_config.h"

//------------------------------------------------------------------------------------------------------

#ifdef CGAL_FOUND

// CGAL needs -DCGAL_NDEBUG to reach peak performance ...
#define CGAL_NDEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Real_timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;

typedef K::Vector_3                               Vector_3;
typedef K::FT                                     FT;
typedef K::Segment_3                              Segment_3;
typedef K::Point_3                                Point_3;

const Point_3 origin = Point_3(CGAL::ORIGIN);

#endif

//------------------------------------------------------------------------------------------------------

#include "atlas/Field.hpp"
#include "atlas/FunctionSpace.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/MeshGen.hpp"
#include "atlas/Parameters.hpp"

namespace atlas {

//------------------------------------------------------------------------------------------------------

const double earthRadius = 6367.47; // from ECMWF model ...

void latlon_to_3d(const double lat, const double lon, double* x, const double r, const double h )
{
    // See http://en.wikipedia.org/wiki/Geodetic_system#From_geodetic_to_ECEF

    double& X = x[0];
    double& Y = x[1];
    double& Z = x[2];

    const double a = r;   // 6378137.0 ;       // WGS84 semi-major axis
    const double e2 = 0;  // ignored -- 6.69437999014E-3; // WGS84 first numerical eccentricity squared

    const double phi = lat / 180.0 * M_PI;
    const double lambda = lon / 180.0 * M_PI;

    const double cos_phi = cos(phi);
    const double sin_phi = sin(phi);
    const double cos_lambda = cos(lambda);
    const double sin_lambda = sin(lambda);

    const double N_phi = a/sqrt(1-e2*sin_phi*sin_phi);

    X = (N_phi + h) * cos_phi * cos_lambda;
    Y = (N_phi + h) * cos_phi * sin_lambda;
    Z = (N_phi * (1-e2) + h) * sin_phi;
}

void latlon_to_3d( const double lat, const double lon, double* x )
{
    latlon_to_3d( lat, lon, x, earthRadius, 0. );
}

//------------------------------------------------------------------------------------------------------

#ifdef CGAL_FOUND


Polyhedron_3* create_convex_hull_from_points( const std::vector< Point3 >& pts )
{
    Polyhedron_3* poly = new Polyhedron_3();

    // insertion from a vector :

    std::vector<Point_3> vertices( pts.size() );
    for( size_t i = 0; i < vertices.size(); ++i )
    {
        vertices[i] = Point_3( pts[i].x[XX], pts[i].x[YY], pts[i].x[ZZ] );
    }

    // compute convex hull of non-collinear points

    CGAL::convex_hull_3( vertices.begin(), vertices.end(), *poly );

    return poly;
}

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

#endif

//------------------------------------------------------------------------------------------------------

atlas::Mesh *atlas::MeshGen::generate_from_ll_points(const std::vector<Point3>& pts)
{
    Mesh* mesh = 0;

#ifdef CGAL_FOUND

    // define polyhedron to hold convex hull

    Polyhedron_3* poly = create_convex_hull_from_points( pts );

    std::cout << "convex hull " << poly->size_of_vertices() << " vertices" << std::endl;

    assert( poly->size_of_vertices() == pts.size() );

    mesh = cgal_polyhedron_to_atlas_mesh( *poly );

    delete poly;
#else

    throw std::string( "CGAL package not found -- triangulation is disabled" );

#endif

    return mesh;
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

