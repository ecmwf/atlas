/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <numeric>

#include "eckit/utils/Hash.h"

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorImpl.h"
#include "atlas/parallel/mpi/mpi.h"

using atlas::Mesh;

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

MeshGeneratorImpl::MeshGeneratorImpl() {}

MeshGeneratorImpl::~MeshGeneratorImpl() {}

Mesh MeshGeneratorImpl::operator()( const Grid& grid ) const {
    Mesh mesh;
    generate( grid, mesh );
    return mesh;
}

Mesh MeshGeneratorImpl::operator()( const Grid& grid, const grid::Distribution& distribution ) const {
    Mesh mesh;
    generate( grid, distribution, mesh );
    return mesh;
}

Mesh MeshGeneratorImpl::generate( const Grid& grid ) const {
    Mesh mesh;
    generate( grid, mesh );
    return mesh;
}

Mesh MeshGeneratorImpl::generate( const Grid& grid, const grid::Distribution& distribution ) const {
    Mesh mesh;
    generate( grid, distribution, mesh );
    return mesh;
}

//----------------------------------------------------------------------------------------------------------------------

void MeshGeneratorImpl::generateGlobalElementNumbering( Mesh& mesh ) const {
    idx_t mpi_size = static_cast<idx_t>( mpi::comm().size() );

    gidx_t loc_nb_elems = mesh.cells().size();
    std::vector<gidx_t> elem_counts( mpi_size );
    std::vector<gidx_t> elem_displs( mpi_size );

    ATLAS_TRACE_MPI( ALLGATHER ) { mpi::comm().allGather( loc_nb_elems, elem_counts.begin(), elem_counts.end() ); }

    elem_displs.at( 0 ) = 0;
    for ( idx_t jpart = 1; jpart < mpi_size; ++jpart ) {
        elem_displs.at( jpart ) = elem_displs.at( jpart - 1 ) + elem_counts.at( jpart - 1 );
    }

    gidx_t gid = 1 + elem_displs.at( mpi::comm().rank() );

    array::ArrayView<gidx_t, 1> glb_idx = array::make_view<gidx_t, 1>( mesh.cells().global_index() );

    for ( idx_t jelem = 0; jelem < mesh.cells().size(); ++jelem ) {
        glb_idx( jelem ) = gid++;
    }

    gidx_t max_glb_idx = std::accumulate( elem_counts.begin(), elem_counts.end(), gidx_t( 0 ) );

    mesh.cells().global_index().metadata().set( "human_readable", true );
    mesh.cells().global_index().metadata().set( "min", 1 );
    mesh.cells().global_index().metadata().set( "max", max_glb_idx );
}

void MeshGeneratorImpl::setProjection( Mesh& mesh, const Projection& p ) const {
    mesh.setProjection( p );
}

void MeshGeneratorImpl::setGrid( Mesh& mesh, const Grid& g, const grid::Distribution& d ) const {
    mesh.setGrid( g );
    mesh.metadata().set( "distribution", d.type() );
}
void MeshGeneratorImpl::setGrid( Mesh& mesh, const Grid& g, const std::string& d ) const {
    mesh.setGrid( g );
    mesh.metadata().set( "distribution", d );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
