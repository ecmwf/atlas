/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>

#include "BuildNode2CellConnectivity.h"

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/domain.h"
#include "atlas/field/Field.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/detail/AccumulateFacets.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/MicroDeg.h"
#include "atlas/util/Unique.h"

using atlas::mesh::detail::accumulate_facets_ordered_by_halo;
using Topology = atlas::mesh::Nodes::Topology;
using atlas::util::microdeg;
using atlas::util::UniqueLonLat;

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

namespace {  // anonymous
struct Sort {
    Sort() = default;
    Sort( gidx_t gid, idx_t idx ) {
        g = gid;
        i = idx;
    }
    gidx_t g;
    idx_t i;
    bool operator<( const Sort& other ) const { return ( g < other.g ); }
};
}  // anonymous namespace

void BuildNode2CellConnectivity::operator()() {
    mesh::Nodes& nodes   = mesh_.nodes();
    const idx_t nb_cells = mesh_.cells().size();

    mesh::Nodes::Connectivity& node_to_cell = nodes.cell_connectivity();
    node_to_cell.clear();

    const mesh::HybridElements::Connectivity& cell_node_connectivity = mesh_.cells().node_connectivity();

    std::vector<idx_t> to_cell_size( nodes.size(), 0 );
    for ( idx_t jcell = 0; jcell < nb_cells; ++jcell ) {
        for ( idx_t j = 0; j < cell_node_connectivity.cols( jcell ); ++j ) {
            ++to_cell_size[cell_node_connectivity( jcell, j )];
        }
    }

    node_to_cell.add( nodes.size(), to_cell_size.data() );

    for ( idx_t jnode = 0; jnode < nodes.size(); ++jnode ) {
        to_cell_size[jnode] = 0;
    }

    UniqueLonLat compute_uid( mesh_ );
    std::vector<Sort> cell_sort( nb_cells );
    for ( idx_t jcell = 0; jcell < nb_cells; ++jcell ) {
        cell_sort[jcell] = Sort( compute_uid( cell_node_connectivity.row( jcell ) ), jcell );
    }


    std::stable_sort( cell_sort.data(), cell_sort.data() + nb_cells );

    for ( idx_t jcell = 0; jcell < nb_cells; ++jcell ) {
        idx_t icell = cell_sort[jcell].i;
        ATLAS_ASSERT( icell < nb_cells );
        for ( idx_t j = 0; j < cell_node_connectivity.cols( icell ); ++j ) {
            idx_t node = cell_node_connectivity( icell, j );
            node_to_cell.set( node, to_cell_size[node]++, icell );
        }
    }
}

BuildNode2CellConnectivity::BuildNode2CellConnectivity( Mesh& mesh ) : mesh_( mesh ) {}


//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {
void atlas__build_node_to_cell_connectivity( Mesh::Implementation* mesh ) {
    ATLAS_ASSERT( mesh != nullptr, "Cannot access uninitialised atlas_Mesh" );
    auto m                               = Mesh( mesh );
    auto build_node_to_cell_connectivity = BuildNode2CellConnectivity( m );
    build_node_to_cell_connectivity();
}
}


//----------------------------------------------------------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
