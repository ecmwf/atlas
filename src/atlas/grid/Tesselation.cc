/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include <vector>
#include <memory>
#include <iostream>

#include "eckit/log/Timer.h"
#include "eckit/memory/ScopedPtr.h"

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

#include "atlas/mesh/Field.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/mesh/FunctionSpace.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Parameters.h"

#include "atlas/grid/PointSet.h"
#include "atlas/grid/Tesselation.h"
#include "atlas/grid/MeshCache.h"
#include "atlas/grid/Grid.h"

using namespace eckit;
using namespace eckit::geometry;

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

#ifdef CGAL_FOUND

Polyhedron_3* create_convex_hull_from_points( const std::vector< Point3 >& pts )
{
    Timer t("convex hull");

    Polyhedron_3* poly = new Polyhedron_3();

    // insertion from a vector :

    std::vector<Point_3> vertices( pts.size() );
    for( size_t i = 0; i < vertices.size(); ++i )
        vertices[i] = Point_3( pts[i](XX), pts[i](YY), pts[i](ZZ) );

    // compute convex hull of non-collinear points

    CGAL::convex_hull_3( vertices.begin(), vertices.end(), *poly );

    return poly;
}

void cgal_polyhedron_to_atlas_mesh(  Mesh& mesh, Polyhedron_3& poly, PointSet& points )
{
    bool ensure_outward_normals = true;

    Timer t ("creating atlas data structure");

    ASSERT( mesh.has_function_space("nodes") );

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT( points.size() == nodes.extents()[0] );

    const size_t nb_nodes = points.size();

    ASSERT( ! mesh.has_function_space("triags") );

    /* triangles */

    const size_t nb_triags = poly.size_of_facets();

    std::vector<int> extents(2);
    extents[0] = nb_triags;
    extents[1] = Field::UNDEF_VARS;

    FunctionSpace& triags  = mesh.add_function_space( new FunctionSpace( "triags", "Lagrange_P1", extents ) );
    triags.metadata().set("type",static_cast<int>(Entity::ELEMS));

    IndexView<int,2> triag_nodes ( triags.create_field<int>("nodes",3) );

    Point3 pt;
    size_t idx[3];
    Polyhedron_3::Vertex_const_handle vts[3];

    std::cout << "inserting triags (" << nb_triags << ")" << std::endl;

//    boost::progress_display show_triag_progress( nb_triags );
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

            pt.assign(p);

            idx[iedge] = points.unique( pt );

            ASSERT( idx[iedge] < nb_nodes );

            vts[iedge] = vh;

            ++iedge;
            ++edge;
        }
        while ( edge != f->facet_begin() && iedge < 3 );

        ASSERT( iedge == 3 );

        if( ensure_outward_normals ) /* ensure outward pointing normal */
        {
            Vector_3 p0 ( origin, vts[0]->point() );
            Vector_3 n  = CGAL::normal( vts[0]->point(), vts[1]->point(), vts[2]->point() );

            FT innerp = n * p0;

            if( innerp < 0 ) // need to swap an edge of the triag
                std::swap( vts[1], vts[2] );
        }

        /* define the triag */

        triag_nodes(tidx,0) = idx[0];
        triag_nodes(tidx,1) = idx[1];
        triag_nodes(tidx,2) = idx[2];

        ++tidx;
//        ++show_triag_progress;
    }

    assert( tidx == nb_triags );
}

#else

struct Polyhedron_3 {};

Polyhedron_3* create_convex_hull_from_points( const std::vector< Point3 >& pts )
{
	throw NotImplemented( "CGAL package not found -- triangulation is disabled", Here() );
}

void cgal_polyhedron_to_atlas_mesh(  Mesh& mesh, Polyhedron_3& poly, PointSet& points )
{
	throw NotImplemented( "CGAL package not found -- triangulation is disabled", Here() );
}

#endif

//------------------------------------------------------------------------------------------------------

void Tesselation::tesselate( Grid& g )
{
	std::string uid = g.uid();

    Mesh& mesh = g.mesh();

	if( MeshCache::get( uid, mesh ) )
        return;

	std::cout << "mesh not in cache -- tesselating grid " << uid << std::endl;

    Tesselation::tesselate( mesh );

	MeshCache::add( uid, mesh );
}

void Tesselation::tesselate( Mesh& mesh )
{
    // don't tesselate meshes already with triags or quads
    if( mesh.has_function_space("triags") || mesh.has_function_space("quads") )
        return;

    Timer t ("grid tesselation");

    // remove duplicate points

    PointSet points( mesh ); /* will remember each point index in readpts */

    std::vector< Point3 > ipts;
    std::vector< size_t > idxs;

    points.list_unique_points( ipts, idxs );

//    std::cout << "unique pts " << ipts.size() << std::endl;
//    std::cout << "duplicates " << points.duplicates().size() << std::endl;

#ifdef CGAL_FOUND

    // define polyhedron to hold convex hull

    eckit::ScopedPtr< Polyhedron_3 > poly( create_convex_hull_from_points( ipts ) );

//    std::cout << "convex hull " << poly->size_of_vertices() << " vertices" << std::endl;

    assert( poly->size_of_vertices() == ipts.size() );

    cgal_polyhedron_to_atlas_mesh( mesh, *poly, points );

#else

	throw NotImplemented( "CGAL package not found -- triangulation is disabled", Here() );

#endif
}

//------------------------------------------------------------------------------------------------------

void Tesselation::create_mesh_structure( Mesh& mesh, const size_t nb_nodes )
{
    // create / ensure mesh has coordinates

    std::vector<int> extents (2);
    if( ! mesh.has_function_space("nodes") )
    {
        extents[0] = nb_nodes;
        extents[1] = Field::UNDEF_VARS;
        FunctionSpace& nodes = mesh.add_function_space( new FunctionSpace( "nodes", "Lagrange_P0", extents ) );
        nodes.metadata().set("type",static_cast<int>(Entity::NODES));
    }

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT(  nodes.extents()[0] == nb_nodes );

    // create / ensure mesh has coordinates

    if( ! nodes.has_field("coordinates") )
        nodes.create_field<double>("coordinates",3);

    // create / ensure mesh has latlon

    if( ! nodes.has_field("latlon") )
        nodes.create_field<double>("latlon",2);

    // create / ensure mesh has global indexes

    if( ! nodes.has_field("glb_idx") )
        nodes.create_field<int>("glb_idx",1);
}

//------------------------------------------------------------------------------------------------------

void Tesselation::generate_latlon_points( Mesh& mesh,
                                          const size_t& nlats,
                                          const size_t& nlong )
{
    const size_t nb_nodes = nlats * nlong;

    Tesselation::create_mesh_structure(mesh,nb_nodes);

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT(  nodes.extents()[0] == nb_nodes );

    ArrayView<double,2> coords  ( nodes.field("coordinates") );
    ArrayView<double,2> latlon  ( nodes.field("latlon") );
    ArrayView<int,   1> glb_idx ( nodes.field("glb_idx") );

    // generate lat/long points

//    std::cout << "generating nlats (" << nlats << ") x  (" << nlong << ")" << " = " << nb_nodes << std::endl;

    const double lat_inc = 180. / (double)nlats;
    const double lat_start = -90 + 0.5*lat_inc;
//    const double lat_end   = 90. - 0.5*lat_inc;

    const double lon_inc = 360. / (double)nlong;
    const double lon_start = 0.5*lon_inc;
//    const double lon_end   = 360. - 0.5*lon_inc;

    double lat = 0;
    double lon = 0;

    size_t visits = 0;

    lat = lat_start;
    for( size_t ilat = 0; ilat < nlats; ++ilat, lat += lat_inc )
    {
        lon = lon_start;
        for( size_t jlon = 0 ; jlon < nlong; ++jlon, lon += lon_inc )
        {
            const size_t idx = jlon + ( ilat * nlong );

            ASSERT( idx < nb_nodes );

            glb_idx(idx) = idx;

            latlon(idx,LAT) = lat;
            latlon(idx,LON) = lon;

            eckit::geometry::latlon_to_3d( lat, lon, coords[idx].data() );

            //            std::cout << idx << " [ " << lat << " ; " << lon << " ] " << p << std::endl;

            ++visits;
        }
    }

    ASSERT( visits == nb_nodes );
}

//------------------------------------------------------------------------------------------------------

void Tesselation::generate_latlon_grid( Mesh& mesh, const size_t& nlats, const size_t& nlong )
{
    const size_t nb_nodes = (nlats+1) * (nlong+1);

    Tesselation::create_mesh_structure(mesh,nb_nodes);

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT( nodes.extents()[0] == nb_nodes );

    ArrayView<double,2> coords  ( nodes.field("coordinates") );
    ArrayView<double,2> latlon  ( nodes.field("latlon") );
    ArrayView<int,   1> glb_idx ( nodes.field("glb_idx") );

    const double lat_inc = 180. / nlats;
    const double lat_start = -90.;
    const double lat_end   =  90.;

    const double lon_inc = 360. / nlong;
    const double lon_start = 0.0;
    const double lon_end   = 360.;

    double lat = 0;
    double lon = 0;

    size_t visits = 0;

    lat = lat_start;
    for( size_t ilat = 0; ilat < nlats+1; ++ilat, lat += lat_inc )
    {
        lon = lon_start;
        for( size_t jlon = 0 ; jlon < nlong+1; ++jlon, lon += lon_inc )
        {

            const size_t idx = jlon + ( ilat * (nlong+1) );

            assert( idx < nb_nodes );

            glb_idx(idx) = idx;

            latlon(idx,LAT) = lat;
            latlon(idx,LON) = lon;

            eckit::geometry::latlon_to_3d( lat, lon, coords[idx].data() );

            ++visits;

//            std::cout << idx << " "
//                      << lat << " "
//                      << lon << " "
//                      << (*pts)[ ilat*(nlats+1) + jlon ](XX) << " "
//                      << (*pts)[ ilat*(nlats+1) + jlon ](YY) << " "
//                      << (*pts)[ ilat*(nlats+1) + jlon ](ZZ) << " "
//                      << std::endl;
//            std::cout << (*pts)[idx] << std::endl;

            if( jlon == nlong ) lon = lon_end;
        }

        if( ilat == nlats ) lat = lat_end;
    }

    ASSERT( visits == nb_nodes );
}

//------------------------------------------------------------------------------------------------------

void Tesselation::create_cell_centres( Mesh& mesh )
{
    ASSERT( mesh.has_function_space("nodes") );
    ASSERT( mesh.has_function_space("triags") );

    FunctionSpace& nodes     = mesh.function_space( "nodes" );
    ArrayView<double,2> coords  ( nodes.field("coordinates") );

    const size_t nb_nodes = nodes.extents()[0];

    FunctionSpace& triags      = mesh.function_space( "triags" );
    IndexView<int,2> triag_nodes ( triags.field( "nodes" ) );

    const size_t nb_triags = triags.extents()[0];

    ArrayView<double,2> triags_centres ( triags.create_field<double>("centre",3) );

    const double third = 1. / 3.;
    for( int e = 0; e < nb_triags; ++e )
    {
        const int i0 =  triag_nodes(e,0);
        const int i1 =  triag_nodes(e,1);
        const int i2 =  triag_nodes(e,2);

        assert( i0 < nb_nodes && i1 < nb_nodes && i2 < nb_nodes );

#if 0 /* print triangle connectivity */
        std::cout << i0 << " " << i1 << " " << i2 << std::endl;
#endif
#if 0 /* print triangle idx and coordinates */
           std::cout << e << " "
                     << i0 << " " << i1 << " " << i2 << " ";
           for( int i = 0; i < 3; ++i )
               std::cout << "("
                     <<  coords(triag_nodes(e,i),XX) << "; "
                     <<  coords(triag_nodes(e,i),YY) << "; "
                     <<  coords(triag_nodes(e,i),ZZ) << ")";
          std::cout << std::endl;
#endif
        triags_centres(e,XX) = third * ( coords(i0,XX) + coords(i1,XX) + coords(i2,XX) );
        triags_centres(e,YY) = third * ( coords(i0,YY) + coords(i1,YY) + coords(i2,YY) );
        triags_centres(e,ZZ) = third * ( coords(i0,ZZ) + coords(i1,ZZ) + coords(i2,ZZ) );

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
        std::cout << triags_centres(e,XX) << " "
                  << triags_centres(e,YY) << " "
                  << triags_centres(e,ZZ) << " "
                  << e << " "
                  << std::endl;
    }
#endif
}

void Tesselation::build_mesh( const Grid& grid, Mesh& mesh )
{
    if( mesh.has_function_space("nodes") ) return;

    const size_t npts = grid.nPoints();

    Tesselation::create_mesh_structure( mesh, npts );

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT(  nodes.extents()[0] == npts );

    ArrayView<double,2> coords  ( nodes.field("coordinates") );
    ArrayView<double,2> latlon  ( nodes.field("latlon") );
    ArrayView<int,   1> glb_idx ( nodes.field("glb_idx") );

    ASSERT( npts == nodes.extents()[0] );

	std::vector<Grid::Point> ll(npts);
	grid.coordinates(ll);

    for( size_t i = 0; i < npts; ++i )
    {
        glb_idx(i) = i;

        double lat = ll[i].lat();
        double lon = ll[i].lon();

        latlon(i,LAT) = lat;
        latlon(i,LON) = lon;

        eckit::geometry::latlon_to_3d( lat, lon, coords[i].data() );
    }
}

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

