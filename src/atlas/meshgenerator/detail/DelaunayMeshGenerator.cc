/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/utils/Hash.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildConvexHull3D.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/mesh/actions/ExtendNodesGlobal.h"
#include "atlas/meshgenerator/detail/DelaunayMeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/projection/Projection.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"

using atlas::Mesh;

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

DelaunayMeshGenerator::DelaunayMeshGenerator() {}

DelaunayMeshGenerator::DelaunayMeshGenerator( const eckit::Parametrisation& ) {}

DelaunayMeshGenerator::~DelaunayMeshGenerator() {}

void DelaunayMeshGenerator::hash( eckit::Hash& h ) const {
    h.add( "Delaunay" );

    // no other settings
}

void DelaunayMeshGenerator::generate( const Grid& grid, const grid::Distribution& dist, Mesh& mesh ) const {
    if ( dist.nb_partitions() > 1 ) {
        Log::warning() << "Delaunay triangulation does not support a GridDistribution"
                          "with more than 1 partition"
                       << std::endl;
        ATLAS_NOTIMPLEMENTED;
        /// TODO: Read mesh on 1 MPI task, and distribute according to
        /// GridDistribution
        /// HINT: use atlas/actions/DistributeMesh
    }
    else {
        generate( grid, mesh );
    }
}

void DelaunayMeshGenerator::generate( const Grid& g, Mesh& mesh ) const {
    createNodes( g, mesh );

    array::ArrayView<gidx_t, 1> gidx = array::make_view<gidx_t, 1>( mesh.nodes().global_index() );
    for ( idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode ) {
        gidx( jnode ) = jnode + 1;
    }

    mesh::actions::BuildXYZField()( mesh );
    mesh::actions::ExtendNodesGlobal()( g,
                                        mesh );  ///< does nothing if global domain
    mesh::actions::BuildConvexHull3D()( mesh );

    setGrid( mesh, g, "serial" );
}

void DelaunayMeshGenerator::createNodes( const Grid& grid, Mesh& mesh ) const {
    idx_t nb_nodes = grid.size();
    mesh.nodes().resize( nb_nodes );

    array::ArrayView<double, 2> xy     = array::make_view<double, 2>( mesh.nodes().xy() );
    array::ArrayView<double, 2> lonlat = array::make_view<double, 2>( mesh.nodes().lonlat() );
    idx_t jnode( 0 );
    Projection projection = grid.projection();
    PointLonLat Pll;
    for ( PointXY Pxy : grid.xy() ) {
        xy( jnode, XX )      = Pxy.x();
        xy( jnode, YY )      = Pxy.y();
        Pll                  = projection.lonlat( Pxy );
        lonlat( jnode, LON ) = Pll.lon();
        lonlat( jnode, LAT ) = Pll.lat();
        ++jnode;
    }
}

namespace {
static MeshGeneratorBuilder<DelaunayMeshGenerator> __delaunay( "delaunay" );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
