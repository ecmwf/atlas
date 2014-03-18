// (C) Copyright 1996-2014 ECMWF.

#include <cmath>
#include <vector>

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

    double& X = x[XX];
    double& Y = x[YY];
    double& Z = x[ZZ];

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
        vertices[i] = Point_3( pts[i](XX), pts[i](YY), pts[i](ZZ) );

//        std::cout << vertices[i] << std::endl;
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

//    std::cout << "nb_nodes = " << nb_nodes << std::endl;

    std::vector<int> bounds(2);
    bounds[0] = Field::UNDEF_VARS;
    bounds[1] = nb_nodes;

    FunctionSpace& nodes = mesh->add_function_space( new FunctionSpace( "nodes", "Lagrange_P0", bounds ) );

    nodes.metadata().set("type",static_cast<int>(Entity::NODES));

    FieldT<double>& coords  = nodes.create_field<double>("coordinates",3);

    FieldT<int>& glb_idx  = nodes.create_field<int>("glb_idx",1);

    std::map< Polyhedron_3::Vertex_const_handle, size_t > vidx;

    size_t inode = 0;
    for( Polyhedron_3::Vertex_const_iterator v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
    {
        vidx[v] = inode;

        glb_idx(inode) = inode;

        const Polyhedron_3::Point_3& p = v->point();

        coords(XX,inode) = p.x();
        coords(YY,inode) = p.y();
        coords(ZZ,inode) = p.z();

//        std::cout << p << std::endl;

        ++inode;
    }

    assert( inode == nb_nodes );

    /* triangles */

    const size_t nb_triags = poly.size_of_facets();

//    std::cout << "nb_triags = " << nb_triags << std::endl;

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

atlas::Mesh* atlas::MeshGen::generate_from_points(const std::vector<Point3>& pts)
{
    Mesh* mesh = 0;

#ifdef CGAL_FOUND

    // define polyhedron to hold convex hull

    Polyhedron_3* poly = create_convex_hull_from_points( pts );

//    std::cout << "convex hull " << poly->size_of_vertices() << " vertices" << std::endl;

    assert( poly->size_of_vertices() == pts.size() );

    mesh = cgal_polyhedron_to_atlas_mesh( *poly );

    delete poly;
#else

    throw std::string( "CGAL package not found -- triangulation is disabled" );

#endif

    return mesh;
}

//------------------------------------------------------------------------------------------------------

std::vector<Point3>* atlas::MeshGen::generate_latlon_points( size_t nlats, size_t nlong )
{
    // generate lat/long points

    const size_t npts = nlats * nlong;

    std::vector< Point3 >* pts = new std::vector< Point3 >( npts );

    const double lat_inc = 180. / nlats;
    const double lat_start = -90 + 0.5*lat_inc;
//    const double lat_end   = 90. - 0.5*lat_inc;

    const double lon_inc = 360. / nlong;
    const double lon_start = 0.5*lon_inc;
//    const double lon_end   = 360. - 0.5*lon_inc;

    double lat = lat_start;
    double lon = lon_start;
    for( size_t ilat = 0; ilat < nlats; ++ilat )
    {
        lon = lon_start;
        for( size_t jlon = 0; jlon < nlong; ++jlon )
        {
            const size_t idx = ilat*nlats + jlon;

            assert( idx < npts );

            atlas::latlon_to_3d( lat, lon, (*pts)[ idx ].data() );

//            std::cout << idx << " "
//                      << lat << " "
//                      << lon << " "
//                      << (*pts)[ ilat*nlats + jlon ].x[XX] << " "
//                      << (*pts)[ ilat*nlats + jlon ].x[YY] << " "
//                      << (*pts)[ ilat*nlats + jlon ].x[ZZ] << " "
//                      << std::endl;
//            std::cout << (*pts)[idx] << std::endl;

            lon += lon_inc;
        }
        lat += lat_inc;
    }

    return pts;
}

//------------------------------------------------------------------------------------------------------

void MeshGen::create_cell_centres(Mesh &mesh)
{
    FunctionSpace& nodes     = mesh.function_space( "nodes" );
    FieldT<double>& coords   = nodes.field<double>( "coordinates" );

    const size_t nb_nodes = nodes.bounds()[1];

    FunctionSpace& triags      = mesh.function_space( "triags" );
    FieldT<int>& triag_nodes   = triags.field<int>( "nodes" );

    const size_t nb_triags = triags.bounds()[1];

    FieldT<double>& triags_centres   = triags.create_field<double>("centre",3);

    const double third = 1. / 3.;
    for( int e = 0; e < nb_triags; ++e )
    {
        const int i0 =  C_IDX( triag_nodes(0,e) );
        const int i1 =  C_IDX( triag_nodes(1,e) );
        const int i2 =  C_IDX( triag_nodes(2,e) );

        assert( i0 < nb_nodes && i1 < nb_nodes && i2 < nb_nodes );

#if 0 /* print triangle connectivity */
        std::cout << i0 << " " << i1 << " " << i2 << std::endl;
#endif
#if 0 /* print triangle idx and coordinates */
           std::cout << e << " "
                     << i0 << " " << i1 << " " << i2 << " ";
           for( int i = 0; i < 3; ++i )
               std::cout << "("
                     <<  coords(XX,C_IDX( triag_nodes(i,e) )) << "; "
                     <<  coords(YY,C_IDX( triag_nodes(i,e) )) << "; "
                     <<  coords(ZZ,C_IDX( triag_nodes(i,e) )) << ")";
          std::cout << std::endl;
#endif
        triags_centres(XX,e) = third * ( coords(XX,i0) + coords(XX,i1) + coords(XX,i2) );
        triags_centres(YY,e) = third * ( coords(YY,i0) + coords(YY,i1) + coords(YY,i2) );
        triags_centres(ZZ,e) = third * ( coords(ZZ,i0) + coords(ZZ,i1) + coords(ZZ,i2) );

#if 0 /* print sorted triangle connectivity */
        std::vector<int> s;
        s.push_back(i0);
        s.push_back(i1);
        s.push_back(i2);
        std::sort(s.begin(),s.end());
        std::cout << s[0] << " " << s[1] << " " << s[2] << std::endl;
#endif
    }

#if 0 /* print triangle baricentres */
    for( int e = 0; e < nb_triags; ++e )
    {
        std::cout << triags_centres(XX,e) << " "
                  << triags_centres(YY,e) << " "
                  << triags_centres(ZZ,e) << " "
                  << e << " "
                  << std::endl;
    }
#endif
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

