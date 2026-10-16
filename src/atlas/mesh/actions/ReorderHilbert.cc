/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/domain/Domain.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/SpaceFillingCurve.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/ReorderHilbert.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace mesh {
namespace actions {

// -------------------------------------------------------------------------------------

// The Hilbert space-filling curve implementation has been moved to a reusable
// component: atlas::grid::HilbertCurve (atlas/grid/SpaceFillingCurve.h), so that
// it can be shared with the space-filling-curve grid partitioner.
using Hilbert = grid::HilbertCurve;

// ------------------------------------------------------------------

ReorderHilbert::ReorderHilbert(const eckit::Parametrisation& config) {
    config.get("recursion", recursion_);
    config.get("ghost_at_end", ghost_at_end_);
}


Domain global_bounding_box(const Mesh& mesh) {
    auto xy = array::make_view<double, 2>(mesh.nodes().xy());

    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();
    for (idx_t i = 0; i < xy.shape(0); ++i) {
        xmin = std::min(xmin, xy(i, XX));
        xmax = std::max(xmax, xy(i, XX));
        ymin = std::min(ymin, xy(i, YY));
        ymax = std::max(ymax, xy(i, YY));
    }
    const auto& comm = atlas::mpi::comm();

    comm.allReduceInPlace(xmin, eckit::mpi::min());
    comm.allReduceInPlace(xmax, eckit::mpi::max());
    comm.allReduceInPlace(ymin, eckit::mpi::min());
    comm.allReduceInPlace(ymax, eckit::mpi::max());
    return RectangularDomain({xmin, xmax}, {ymin, ymax});
}

std::vector<idx_t> ReorderHilbert::computeNodesOrder(Mesh& mesh) {
    using hilbert_reordering_t = std::vector<std::pair<gidx_t, idx_t>>;

    Hilbert hilbert{global_bounding_box(mesh), recursion_};

    auto xy    = array::make_view<double, 2>(mesh.nodes().xy());
    auto ghost = array::make_view<int, 1>(mesh.nodes().ghost());

    idx_t size = xy.shape(0);
    hilbert_reordering_t hilbert_reordering;
    hilbert_reordering.reserve(size);
    ATLAS_TRACE_SCOPE("hilbert nodes") {
        for (idx_t n = 0; n < size; ++n) {
            PointXY p{xy(n, XX), xy(n, YY)};
            if (not ghost(n)) {
                hilbert_reordering.emplace_back(hilbert(p), n);
            }
            else {
                if (ghost_at_end_) {
                    // ghost nodes get a fake "hilbert_idx" at the end
                    hilbert_reordering.emplace_back(hilbert.nb_keys() + n, n);
                }
                else {
                    hilbert_reordering.emplace_back(hilbert(p), n);
                }
            }
        }
    }

    std::sort(hilbert_reordering.begin(), hilbert_reordering.end());
    std::vector<idx_t> order;
    order.reserve(size);
    for (const auto& pair : hilbert_reordering) {
        order.emplace_back(pair.second);
    }
    return order;
}


#if 0
// Reorder elements
if ( 0 ) {
    hilbert_reordering_t hilbert_reordering;
    std::vector<idx_t> order;
    std::vector<idx_t> order_inverse;

    auto cell_centres =
        Field( "cell_centres", array::make_datatype<double>(), array::make_shape( mesh.cells().size(), 2 ) );
    auto nodes_xy = array::make_view<double, 2>( mesh.nodes().xy() );
    for ( idx_t t = 0; t < mesh.cells().nb_types(); ++t ) {
        auto& cells = mesh.cells().elements( t );
        auto xy     = cells.view<double, 2>( cell_centres );
        auto flags  = cells.view<int, 1>( mesh.cells().flags() );
        auto halo   = cells.view<idx_t, 1>( mesh.cells().halo() );

        // Compute cell-centres
        {
            const auto& node_connectivity = cells.node_connectivity();
            const idx_t nb_nodes          = cells.nb_nodes();
            const double nb_nodes_double  = nb_nodes;
            for ( idx_t e = 0; e < cells.size(); ++e ) {
                double x{0};
                double y{0};
                for ( idx_t c = 0; c < nb_nodes; ++c ) {
                    idx_t n = node_connectivity( e, c );
                    x += nodes_xy( n, XX );
                    y += nodes_xy( n, YY );
                }
                xy( e, XX ) = x / nb_nodes_double;
                xy( e, YY ) = y / nb_nodes_double;
            }
        }


        auto skip = [&]( idx_t n ) {
            if ( halo( n ) || mesh::Nodes::Topology::check( flags( n ), mesh::Nodes::Topology::PATCH ) ) {
                return true;
            }
            return false;
        };
        idx_t size = xy.shape( 0 );
        hilbert_reordering.clear();
        hilbert_reordering.reserve( size );
        ATLAS_TRACE_SCOPE( "hilbert elements[" + std::to_string( t ) + "]" ) {
            for ( idx_t n = 0; n < size; ++n ) {
                PointXY p{xy( n, XX ), xy( n, YY )};
                if ( not skip( n ) ) {
                    hilbert_reordering.emplace_back( hilbert( p ), n );
                }
                else {  // halo elements get a fake "hilbert_idx" at the end
                    hilbert_reordering.emplace_back( hilbert.nb_keys() + n, n );
                }
            }
        }
        ATLAS_ASSERT( hilbert_reordering.size() == size );
        std::sort( hilbert_reordering.begin(), hilbert_reordering.end() );
        order.clear();
        order.reserve( size );
        order_inverse.resize( size );
        idx_t c{0};
        for ( const auto& pair : hilbert_reordering ) {
            order.emplace_back( pair.second );
            order_inverse[pair.second] = c++;
            ATLAS_ASSERT( pair.second < size );
        }

        for ( idx_t ifield = 0; ifield < mesh.cells().nb_fields(); ++ifield ) {
            reorder_field( mesh.cells().field( ifield ), order, cells.begin(), cells.end() );
        }

        reorder_connectivity( cells.node_connectivity(), order );
    }
}
#endif

namespace {
static ReorderBuilder<ReorderHilbert> __ReorderHilbert("hilbert");
}  // namespace


}  // namespace actions
}  // namespace mesh
}  // namespace atlas
