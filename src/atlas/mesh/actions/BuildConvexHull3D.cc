/*
 * (C) Copyright 1996-2017 ECMWF.
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
#include "eckit/log/BigNum.h"
#include "eckit/memory/ScopedPtr.h"

#include "atlas/library/config.h"

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

#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/actions/BuildConvexHull3D.h"
#include "atlas/field/Field.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/interpolation/method/PointSet.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"

using namespace eckit;
using namespace eckit::geometry;

using atlas::interpolation::method::PointSet;
using atlas::interpolation::method::PointIndex3;

namespace atlas {
namespace mesh {
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

    mesh::Nodes& nodes = mesh.nodes();

    ASSERT( points.size() == nodes.size() );

    const idx_t nb_nodes = idx_t(points.size());

    ASSERT( mesh.cells().size() == 0 );

    /* triangles */

    const size_t nb_triags = poly.size_of_facets();
    mesh.cells().add( new mesh::temporary::Triangle(), nb_triags );
    mesh::HybridElements::Connectivity& triag_nodes = mesh.cells().node_connectivity();
    array::ArrayView<gidx_t,1> triag_gidx = array::make_view<gidx_t,1>( mesh.cells().global_index() );
    array::ArrayView<int,1> triag_part    = array::make_view<int,1>( mesh.cells().partition() );

    Point3 pt;
    idx_t idx[3];
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

        triag_nodes.set(tidx, idx );

        triag_gidx(tidx) = tidx+1;

        triag_part(tidx) = 0;

        ++tidx;
    }

    ASSERT( tidx == nb_triags );
}

#else

struct Polyhedron_3 {
  size_t size_of_vertices() const { return 0; }
};

static Polyhedron_3* create_convex_hull_from_points( const std::vector< Point3 >& pts )
{
  throw NotImplemented( "CGAL package not found -- Delaunay triangulation is disabled", Here() );
}

static void cgal_polyhedron_to_atlas_mesh(  Mesh& mesh, Polyhedron_3& poly, PointSet& points )
{
  throw NotImplemented( "CGAL package not found -- Delaunay triangulation is disabled", Here() );
}

#endif

//----------------------------------------------------------------------------------------------------------------------

void BuildConvexHull3D::operator()( Mesh& mesh ) const
{
    // don't tesselate meshes already with triags or quads
    if( mesh.cells().size() )
        return;

    Timer t ("grid tesselation");

    // remove duplicate points

    PointSet points( mesh );

    std::vector< Point3 > ipts;

    points.list_unique_points( ipts );

//    std::cout << "unique pts " << ipts.size() << std::endl;
//    std::cout << "duplicates " << points.duplicates().size() << std::endl;


    // define polyhedron to hold convex hull

    eckit::ScopedPtr< Polyhedron_3 > poly( create_convex_hull_from_points( ipts ) );

//    std::cout << "convex hull " << poly->size_of_vertices() << " vertices" << std::endl;

    ASSERT( poly->size_of_vertices() == ipts.size() );

    cgal_polyhedron_to_atlas_mesh( mesh, *poly, points );

}

//----------------------------------------------------------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas


