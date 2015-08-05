/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/actions/BuildConvexHull3D.h"

#include <cmath>
#include <vector>
#include <memory>
#include <iostream>

#include "eckit/config/Resource.h"
#include "eckit/log/Timer.h"
#include "eckit/log/BigNum.h"

#include "eckit/memory/ScopedPtr.h"

#include "atlas/atlas_config.h"

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

#include "atlas/Field.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/Nodes.h"
#include "atlas/Parameters.h"
#include "atlas/PointSet.h"
#include "atlas/Grid.h"

using namespace eckit;
using namespace eckit::geometry;

namespace atlas {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

#ifdef CGAL_FOUND

static Polyhedron_3* create_convex_hull_from_points( const std::vector< Point3 >& pts )
{
    Timer t("Convex hull");

    Polyhedron_3* poly = new Polyhedron_3();

    // insertion from a vector :

    std::vector<Point_3> vertices( pts.size() );
    for( size_t i = 0; i < vertices.size(); ++i )
        vertices[i] = Point_3( pts[i](XX), pts[i](YY), pts[i](ZZ) );

    // compute convex hull of non-collinear points

    CGAL::convex_hull_3( vertices.begin(), vertices.end(), *poly );

    return poly;
}

static void cgal_polyhedron_to_atlas_mesh(  Mesh& mesh, Polyhedron_3& poly, PointSet& points )
{
    bool ensure_outward_normals = true;

    Timer t ("Creating atlas data structure");

    Nodes& nodes = mesh.nodes();

    ASSERT( points.size() == nodes.size() );

    const size_t nb_nodes = points.size();

    ASSERT( ! mesh.has_function_space("triags") );

    /* triangles */

    const size_t nb_triags = poly.size_of_facets();

    std::vector<size_t> extents(2);
    extents[0] = nb_triags;
    extents[1] = FunctionSpace::UNDEF_VARS;

    FunctionSpace& triags  = mesh.create_function_space( "triags", "Lagrange_P1", extents );
    triags.metadata().set<long>("type",Entity::ELEMS);

    IndexView<int,2> triag_nodes   ( triags.create_field<int>("nodes",3) );
    ArrayView<gidx_t,1> triag_gidx ( triags.create_field<gidx_t>("glb_idx",1) );
    ArrayView<int,1> triag_part    ( triags.create_field<int>("partition",1) );

    Point3 pt;
    size_t idx[3];
    Polyhedron_3::Vertex_const_handle vts[3];

    std::cout << "Inserting triags (" << eckit::BigNum(nb_triags) << ")" << std::endl;

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

        triag_gidx(tidx) = tidx+1;

        triag_part(tidx) = 0;

        ++tidx;
    }

    ASSERT( tidx == nb_triags );

    /* quads */

    const size_t nb_quads = 0;

    std::vector<size_t> quads_extents(2);
    extents[0] = nb_quads;
    extents[1] = FunctionSpace::UNDEF_VARS;

    FunctionSpace& quads  = mesh.create_function_space( "quads", "Lagrange_P1", quads_extents );
    quads.metadata().set<long>("type",Entity::ELEMS);

    IndexView<int,2> quads_nodes   ( quads.create_field<int>("nodes",3) );
    ArrayView<gidx_t,1> quads_gidx ( quads.create_field<gidx_t>("glb_idx",1) );
    ArrayView<int,1> quads_part    ( quads.create_field<int>("partition",1) );
}

#else

struct Polyhedron_3 {};

static Polyhedron_3* create_convex_hull_from_points( const std::vector< Point3 >& pts )
{
	throw NotImplemented( "CGAL package not found -- triangulation is disabled", Here() );
}

static void cgal_polyhedron_to_atlas_mesh(  Mesh& mesh, Polyhedron_3& poly, PointSet& points )
{
	throw NotImplemented( "CGAL package not found -- triangulation is disabled", Here() );
}

#endif

//----------------------------------------------------------------------------------------------------------------------

void BuildConvexHull3D::operator()( Mesh& mesh ) const
{
    // don't tesselate meshes already with triags or quads
    if( mesh.has_function_space("triags") || mesh.has_function_space("quads") )
        return;

    Timer t ("grid tesselation");

    // remove duplicate points

    PointSet points( mesh );

    std::vector< Point3 > ipts;

    points.list_unique_points( ipts );

//    std::cout << "unique pts " << ipts.size() << std::endl;
//    std::cout << "duplicates " << points.duplicates().size() << std::endl;

#ifdef CGAL_FOUND

    // define polyhedron to hold convex hull

    eckit::ScopedPtr< Polyhedron_3 > poly( create_convex_hull_from_points( ipts ) );

//    std::cout << "convex hull " << poly->size_of_vertices() << " vertices" << std::endl;

    ASSERT( poly->size_of_vertices() == ipts.size() );

    cgal_polyhedron_to_atlas_mesh( mesh, *poly, points );

#else

        throw NotImplemented( "CGAL package not found -- Delaunay triangulation is disabled", Here() );

#endif
}

//----------------------------------------------------------------------------------------------------------------------

} // actions
} // atlas


