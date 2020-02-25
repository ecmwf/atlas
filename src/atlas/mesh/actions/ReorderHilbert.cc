/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <array>
#include <bitset>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/domain/Domain.h"
#include "atlas/grid/Grid.h"
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

/// @brief Class to compute a global index given a coordinate, based on the
/// Hilbert Spacefilling Curve.
///
/// This algorithm is based on:
/// - John J. Bartholdi and Paul Goldsman "Vertex-Labeling Algorithms for the Hilbert Spacefilling Curve"\n
/// It is adapted to return contiguous numbers of the gidx_t type, instead of a double [0,1]
///
/// Given a bounding box and number of hilbert recursions, the bounding box can be divided in
/// 2^(dim*levels) equally spaced cells. A given coordinate falling inside one of these cells, is assigned
/// the 1-dimensional Hilbert-index of this cell. To make sure that 1 coordinate corresponds to only 1
/// Hilbert index, the number of levels have to be increased.
/// In 2D, the recursion cannot be higher than 15, if you want the indices to fit in "unsigned int" type of 32bit.
/// In 2D, the recursion cannot be higher than 30, if you want the indices to fit in "unsigned int" type of 64bit.
///
///
/// No attempt is made to provide the most efficient algorithm. There exist other open-source
/// libraries with more efficient algorithms, such as libhilbert, but its LGPL license
/// is not compatible with this licence.
///
/// @author Willem Deconinck
class Hilbert {
public:
    /// Constructor
    /// Initializes the hilbert space filling curve with a given "space" and "levels"
    Hilbert( const Domain& domain, idx_t levels );

    /// Compute the hilbert code for a given point in 2D
    gidx_t operator()( const PointXY& point );

    /// Compute the hilbert code for a given point in 2D
    /// @param [out] relative_tolerance  cell-size of smallest level divided by bounding-box size
    gidx_t operator()( const PointXY& point, double& relative_tolerance );

    /// Return the maximum hilbert code possible with the initialized levels
    ///
    /// Care has to be taken that this number is not larger than the precision of the type storing
    /// the hilbert codes.
    gidx_t nb_keys() const { return nb_keys_; }

private:  // functions
    using box_t = std::array<PointXY, 4>;

    /// @brief Recursive algorithm
    gidx_t recursive_algorithm( const PointXY& p, const box_t& box, idx_t level );

private:  // data
    /// Vertex label type (4 vertices in 2D)
    enum VertexLabel
    {
        A = 0,
        B = 1,
        C = 2,
        D = 3
    };

    /// Bounding box, defining the space to be filled
    const RectangularDomain domain_;

    /// maximum recursion level of the Hilbert space filling curve
    idx_t max_level_;

    /// maximum number of unique codes, computed by max_level
    gidx_t nb_keys_;
    gidx_t nb_keys_2_;
};

// -------------------------------------------------------------------------------------

Hilbert::Hilbert( const Domain& domain, idx_t levels ) : domain_{domain}, max_level_( levels ) {
    nb_keys_2_ = gidx_t( std::pow( gidx_t( 4 ), gidx_t( max_level_ ) ) );
    nb_keys_   = nb_keys_2_ * 2;
}


gidx_t Hilbert::operator()( const PointXY& point ) {
    box_t box;
    box[A]            = {domain_.xmin(), domain_.ymax()};
    box[B]            = {domain_.xmin(), domain_.ymin()};
    box[C]            = {domain_.xmax(), domain_.ymin()};
    box[D]            = {domain_.xmax(), domain_.ymax()};
    const double xmid = ( domain_.xmin() + domain_.xmax() ) * 0.5;
    if ( point.x() < xmid ) {
        box[C].x() = xmid;
        box[D].x() = xmid;
        return recursive_algorithm( point, box, 0 );
    }
    else {
        box[A].x() = xmid;
        box[B].x() = xmid;
        return recursive_algorithm( point, box, 0 ) + nb_keys_2_;
    }
}

gidx_t Hilbert::recursive_algorithm( const PointXY& p, const box_t& box, idx_t level ) {
    if ( level == max_level_ ) {
        return 0;
    }

    double min_distance = std::numeric_limits<double>::max();

    auto compute_distance2 = []( const PointXY& p1, const PointXY& p2 ) {
        // workaround because of eckit 1.3.2 issue with constness in KPoint
        double d = 0;
        for ( size_t i = 0; i < 2; i++ ) {
            double dx = p1[i] - p2[i];
            d += dx * dx;
        }
        return d;
    };

    auto compute_average = []( const PointXY& p1, const PointXY& p2 ) {
        // workaround because of eckit 1.3.2 issue with constness in KPoint
        PointXY avg;
        avg.x() = p1.x() + p2.x();
        avg.x() *= 0.5;
        avg.y() = p1.y() + p2.y();
        avg.y() *= 0.5;
        return avg;
    };

    idx_t quadrant{0};
    for ( idx_t idx = 0; idx < 4; ++idx ) {
        // double distance = box[idx].distance2( p );  // does not compile with eckit 1.3.2
        double distance = compute_distance2( p, box[idx] );  // workaround
        if ( distance < min_distance ) {
            quadrant     = idx;
            min_distance = distance;
        }
    }

    box_t box_quadrant;
    switch ( quadrant ) {
        case A:
            box_quadrant[A] = box[A];
            // box_quadrant[B] = ( box[A] + box[D] ) * 0.5;  // does not compile with eckit 1.3.2
            // box_quadrant[C] = ( box[A] + box[C] ) * 0.5;  // does not compile with eckit 1.3.2
            // box_quadrant[D] = ( box[A] + box[B] ) * 0.5;  // does not compile with eckit 1.3.2
            box_quadrant[B] = compute_average( box[A], box[D] );  // workaround
            box_quadrant[C] = compute_average( box[A], box[C] );  // workaround
            box_quadrant[D] = compute_average( box[A], box[B] );  // workaround
            break;
        case B:
            // box_quadrant[A] = ( box[B] + box[A] ) * 0.5;  // does not compile with eckit 1.3.2
            box_quadrant[B] = box[B];
            // box_quadrant[C] = ( box[B] + box[C] ) * 0.5;  // does not compile with eckit 1.3.2
            // box_quadrant[D] = ( box[B] + box[D] ) * 0.5;  // does not compile with eckit 1.3.2
            box_quadrant[A] = compute_average( box[B], box[A] );  // workaround
            box_quadrant[C] = compute_average( box[B], box[C] );  // workaround
            box_quadrant[D] = compute_average( box[B], box[D] );  // workaround
            break;
        case C:
            // box_quadrant[A] = ( box[C] + box[A] ) * 0.5;  // does not compile with eckit 1.3.2
            // box_quadrant[B] = ( box[C] + box[B] ) * 0.5;  // does not compile with eckit 1.3.2
            box_quadrant[C] = box[C];
            // box_quadrant[D] = ( box[C] + box[D] ) * 0.5;  // does not compile with eckit 1.3.2
            box_quadrant[A] = compute_average( box[C], box[A] );  // workaround
            box_quadrant[B] = compute_average( box[C], box[B] );  // workaround
            box_quadrant[D] = compute_average( box[C], box[D] );  // workaround

            break;
        case D:
            // box_quadrant[A] = ( box[D] + box[C] ) * 0.5;  // does not compile with eckit 1.3.2
            // box_quadrant[B] = ( box[D] + box[B] ) * 0.5;  // does not compile with eckit 1.3.2
            // box_quadrant[C] = ( box[D] + box[A] ) * 0.5;  // does not compile with eckit 1.3.2
            box_quadrant[D] = box[D];
            box_quadrant[A] = compute_average( box[D], box[C] );  // workaround
            box_quadrant[B] = compute_average( box[D], box[B] );  // workaround
            box_quadrant[C] = compute_average( box[D], box[A] );  // workaround

            break;
    }

    // The key has 4 possible values per recursion (1 for each quadrant),
    // which can be represented by 2 bits per recursion
    //   A --> 00
    //   B --> 01
    //   C --> 10
    //   D --> 11
    // Trailing zero-bits are added depending on the level:
    //   level max_level_-1 --> none
    //   level max_level_-2 --> 00
    //   level max_level_-2 --> 0000
    //   level max_level_-3 --> 000000
    gidx_t key = 0;
    auto index = ( max_level_ - level ) * 2 - 1;
    gidx_t mask;

    // Create a mask value with all trailing bits for leftmost bit (of 2)
    mask = gidx_t( 1 ) << index;

    // Add mask to key
    if ( quadrant == C || quadrant == D ) {
        key |= mask;
    }

    // Create a mask value with all trailing bits for rightmost bit (of 2)
    mask = gidx_t( 1 ) << ( index - 1 );

    // Add mask to key
    if ( quadrant == B || quadrant == D ) {
        key |= mask;
    }

    return recursive_algorithm( p, box_quadrant, level + 1 ) + key;
}

// ------------------------------------------------------------------

ReorderHilbert::ReorderHilbert( const eckit::Parametrisation& config ) {
    config.get( "recursion", recursion_ );
    config.get( "ghost_at_end", ghost_at_end_ );
}


Domain global_bounding_box( const Mesh& mesh ) {
    auto xy = array::make_view<double, 2>( mesh.nodes().xy() );

    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();
    for ( idx_t i = 0; i < xy.shape( 0 ); ++i ) {
        xmin = std::min( xmin, xy( i, XX ) );
        xmax = std::max( xmax, xy( i, XX ) );
        ymin = std::min( ymin, xy( i, YY ) );
        ymax = std::max( ymax, xy( i, YY ) );
    }
    const auto& comm = atlas::mpi::comm();

    comm.allReduceInPlace( xmin, eckit::mpi::min() );
    comm.allReduceInPlace( xmax, eckit::mpi::max() );
    comm.allReduceInPlace( ymin, eckit::mpi::min() );
    comm.allReduceInPlace( ymax, eckit::mpi::max() );
    return RectangularDomain( {xmin, xmax}, {ymin, ymax} );
}

std::vector<idx_t> ReorderHilbert::computeNodesOrder( Mesh& mesh ) {
    using hilbert_reordering_t = std::vector<std::pair<gidx_t, idx_t>>;

    Hilbert hilbert{global_bounding_box( mesh ), recursion_};

    auto xy    = array::make_view<double, 2>( mesh.nodes().xy() );
    auto ghost = array::make_view<int, 1>( mesh.nodes().ghost() );

    idx_t size = xy.shape( 0 );
    hilbert_reordering_t hilbert_reordering;
    hilbert_reordering.reserve( size );
    ATLAS_TRACE_SCOPE( "hilbert nodes" ) {
        for ( idx_t n = 0; n < size; ++n ) {
            PointXY p{xy( n, XX ), xy( n, YY )};
            if ( not ghost( n ) ) {
                hilbert_reordering.emplace_back( hilbert( p ), n );
            }
            else {
                if ( ghost_at_end_ ) {
                    // ghost nodes get a fake "hilbert_idx" at the end
                    hilbert_reordering.emplace_back( hilbert.nb_keys() + n, n );
                }
                else {
                    hilbert_reordering.emplace_back( hilbert( p ), n );
                }
            }
        }
    }

    std::sort( hilbert_reordering.begin(), hilbert_reordering.end() );
    std::vector<idx_t> order;
    order.reserve( size );
    for ( const auto& pair : hilbert_reordering ) {
        order.emplace_back( pair.second );
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
static ReorderBuilder<ReorderHilbert> __ReorderHilbert( "hilbert" );
}  // namespace


}  // namespace actions
}  // namespace mesh
}  // namespace atlas
