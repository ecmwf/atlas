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

#include "eckit/config/Resource.h"
#include "eckit/log/Timer.h"
#include "eckit/log/BigNum.h"

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

#include "atlas/Field.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/Parameters.h"

#include "atlas/PointSet.h"
#include "atlas/Tesselation.h"
#include "atlas/MeshCache.h"
#include "atlas/Grid.h"

#include "atlas/grids/ReducedGrid.h"

#include "atlas/meshgen/ReducedGridMeshGenerator.h"

using namespace eckit;
using namespace eckit::geometry;

namespace atlas {


//------------------------------------------------------------------------------------------------------

#ifdef CGAL_FOUND

Polyhedron_3* create_convex_hull_from_points( const std::vector< Point3 >& pts )
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

void cgal_polyhedron_to_atlas_mesh(  Mesh& mesh, Polyhedron_3& poly, PointSet& points )
{
    bool ensure_outward_normals = true;

    Timer t ("Creating atlas data structure");

    ASSERT( mesh.has_function_space("nodes") );

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT( points.size() == nodes.shape(0) );

    const size_t nb_nodes = points.size();

    ASSERT( ! mesh.has_function_space("triags") );

    /* triangles */

    const size_t nb_triags = poly.size_of_facets();

    std::vector<int> extents(2);
    extents[0] = nb_triags;
    extents[1] = Field::UNDEF_VARS;

    FunctionSpace& triags  = mesh.create_function_space( "triags", "Lagrange_P1", extents );
    triags.metadata().set("type",static_cast<int>(Entity::ELEMS));

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

    std::vector<int> quads_extents(2);
    extents[0] = nb_quads;
    extents[1] = Field::UNDEF_VARS;

    FunctionSpace& quads  = mesh.create_function_space( "quads", "Lagrange_P1", quads_extents );
    quads.metadata().set("type",static_cast<int>(Entity::ELEMS));

    IndexView<int,2> quads_nodes   ( quads.create_field<int>("nodes",3) );
    ArrayView<gidx_t,1> quads_gidx ( quads.create_field<gidx_t>("glb_idx",1) );
    ArrayView<int,1> quads_part    ( quads.create_field<int>("partition",1) );
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

/// @ TODO Abstract this method into a MeshGenerator class (with use of MeshCache)

void Tesselation::tesselate(const Grid& g, Mesh& mesh) {

  std::string uid = g.unique_id();

  MeshCache cache;

  if (cache.retrieve(g, mesh)) return;

  std::cout << "Mesh not in cache -- tesselating grid " << uid << std::endl;

  bool atlasTriangulateRG = eckit::Resource<bool>("$ATLAS_TRIANGULATE_RG", false);

  const grids::ReducedGrid* rg = dynamic_cast<const grids::ReducedGrid*>(&g);
  if (atlasTriangulateRG && rg) {

    // fast tesselation method, specific for ReducedGrid's

    std::cout << "Mesh is ReducedGrid " << g.shortName() << std::endl;

    ASSERT(rg);

    meshgen::ReducedGridMeshGenerator mg;

    // force these flags
    mg.options.set("three_dimensional",true);
    mg.options.set("patch_pole",true);
    mg.options.set("include_pole",false);
    mg.options.set("triangulate",true);

    mg.generate(*rg, mesh);

  } else {

    // slower, more robust tesselation method, using Delaunay triangulation

    std::cout << "Using Delaunay triangulation on grid: " << g.shortName() << std::endl;

    Tesselation::delaunay_triangulation(mesh);
  }

  cache.insert(g, mesh);

}

void Tesselation::delaunay_triangulation( Mesh& mesh )
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

        throw NotImplemented( "CGAL package not found -- Delaunay triangulation is disabled", Here() );

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
		FunctionSpace& nodes = mesh.create_function_space(  "nodes", "Lagrange_P0", extents );
        nodes.metadata().set("type",static_cast<int>(Entity::NODES));
    }

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT(  nodes.shape(0) == nb_nodes );

    // create / ensure mesh has coordinates

	nodes.create_field<double>("xyz",3,IF_EXISTS_RETURN);

    // create / ensure mesh has latlon

	nodes.create_field<double>("lonlat",2,IF_EXISTS_RETURN);

    // create / ensure mesh has global indexes

	nodes.create_field<gidx_t>("glb_idx",1,IF_EXISTS_RETURN);
}

//------------------------------------------------------------------------------------------------------

void Tesselation::create_cell_centres( Mesh& mesh )
{
    ASSERT( mesh.has_function_space("nodes") );
    ASSERT( mesh.has_function_space("triags") );
    ASSERT( mesh.has_function_space("quads") );

    FunctionSpace& nodes     = mesh.function_space( "nodes" );
    ArrayView<double,2> coords  ( nodes.field("xyz") );

    const size_t nb_nodes = nodes.shape(0);

    if( mesh.has_function_space("triags") ) {

        FunctionSpace& triags      = mesh.function_space( "triags" );
        IndexView<int,2> triag_nodes ( triags.field( "nodes" ) );
        const size_t nb_triags = triags.shape(0);

        ArrayView<double,2> triags_centres ( triags.create_field<double>("centre",3) );

        const double third = 1. / 3.;
        for( int e = 0; e < nb_triags; ++e )
        {
            const int i0 =  triag_nodes(e,0);
            const int i1 =  triag_nodes(e,1);
            const int i2 =  triag_nodes(e,2);

            assert( i0 < nb_nodes && i1 < nb_nodes && i2 < nb_nodes );

            triags_centres(e,XX) = third * ( coords(i0,XX) + coords(i1,XX) + coords(i2,XX) );
            triags_centres(e,YY) = third * ( coords(i0,YY) + coords(i1,YY) + coords(i2,YY) );
            triags_centres(e,ZZ) = third * ( coords(i0,ZZ) + coords(i1,ZZ) + coords(i2,ZZ) );

        }
    }

    if( mesh.has_function_space("quads") ) {
        FunctionSpace& quads  = mesh.function_space( "quads" );
        IndexView<int,2> quads_nodes ( quads.field( "nodes" ) );
        const size_t nb_quads = quads.shape(0);

        ArrayView<double,2> quads_centres ( quads.create_field<double>("centre",3) );

        const double fourth = 1. / 4.;
        for( int e = 0; e < nb_quads; ++e )
        {
            const int i0 =  quads_nodes(e,0);
            const int i1 =  quads_nodes(e,1);
            const int i2 =  quads_nodes(e,2);
            const int i3 =  quads_nodes(e,3);

            assert( i0 < nb_nodes && i1 < nb_nodes && i2 < nb_nodes && i3 < nb_nodes );

            quads_centres(e,XX) = fourth * ( coords(i0,XX) + coords(i1,XX) + coords(i2,XX) + coords(i3,XX) );
            quads_centres(e,YY) = fourth * ( coords(i0,YY) + coords(i1,YY) + coords(i2,YY) + coords(i3,YY) );
            quads_centres(e,ZZ) = fourth * ( coords(i0,ZZ) + coords(i1,ZZ) + coords(i2,ZZ) + coords(i3,ZZ) );

        }
    }
}

void Tesselation::build_mesh( const Grid& grid, Mesh& mesh )
{
    if( mesh.has_function_space("nodes") ) return;

    const size_t npts = grid.npts();

    Tesselation::create_mesh_structure( mesh, npts );

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT(  nodes.shape(0) == npts );

    ArrayView<double,2> coords  ( nodes.field("xyz") );
    ArrayView<double,2> lonlat  ( nodes.field("lonlat") );
    ArrayView<gidx_t,1> glb_idx ( nodes.field("glb_idx") );

    ASSERT( npts == nodes.shape(0) );

    std::vector<Grid::Point> ll(npts);
    grid.lonlat(lonlat.data());

    for( size_t i = 0; i < npts; ++i )
    {
        glb_idx(i) = i;

        double lon = lonlat(i,LON);
        double lat = lonlat(i,LAT);

        eckit::geometry::lonlat_to_3d( lon, lat, coords[i].data() );
    }
}

//------------------------------------------------------------------------------------------------------


} // namespace atlas

