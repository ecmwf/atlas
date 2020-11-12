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

#include "atlas/mesh/actions/Reorder.h"

// For static linking
#include "atlas/mesh/actions/ReorderHilbert.h"
#include "atlas/mesh/actions/ReorderReverseCuthillMckee.h"


#include "atlas/array.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

namespace atlas {
namespace mesh {
namespace actions {

// ------------------------------------------------------------------

struct force_link {
    template <typename T>
    void load_builder() {
        ReorderBuilder<T>( "tmp" );
    }
    force_link() {
        load_builder<ReorderHilbert>();
        load_builder<ReorderReverseCuthillMckee>();
    }
};


// ------------------------------------------------------------------

const ReorderImpl* ReorderFactory::build( const eckit::Parametrisation& config ) {
    static force_link static_linking;
    std::string builder{"none"};
    config.get( "type", builder );
    auto factory = get( builder );
    return factory->make( config );
}

// ------------------------------------------------------------------

template <typename Value, int Rank>
struct ReorderField {};

// ------------------------------------------------------------------

template <typename Value>
struct ReorderField<Value, 1> {
    static constexpr int Rank = 1;
    static std::string apply( Field& field, const std::vector<idx_t>& order, idx_t begin, idx_t end ) {
        auto array = array::make_view<Value, Rank>( field );
        end        = std::min( end, array.shape( 0 ) );
        idx_t size = end - begin;
        array::ArrayT<Value> tmp_array( size );
        auto tmp = array::make_view<Value, Rank>( tmp_array );
        for ( idx_t n = 0; n < size; ++n ) {
            tmp( n ) = array( begin + n );
        }
        for ( idx_t n = 0; n < size; ++n ) {
            array( begin + n ) = tmp( order[n] );
        }
        return field.name();
    }
};

// ------------------------------------------------------------------

template <typename Value>
struct ReorderField<Value, 2> {
    static constexpr int Rank = 2;
    static std::string apply( Field& field, const std::vector<idx_t>& order, idx_t begin, idx_t end ) {
        auto array = array::make_view<Value, Rank>( field );
        end        = std::min( end, array.shape( 0 ) );
        idx_t size = end - begin;
        array::ArrayT<Value> tmp_array( size, field.shape( 1 ) );
        auto tmp = array::make_view<Value, Rank>( tmp_array );
        for ( idx_t n = 0; n < size; ++n ) {
            for ( idx_t v = 0; v < array.shape( 1 ); ++v ) {
                tmp( n, v ) = array( begin + n, v );
            }
        }
        for ( idx_t n = 0; n < size; ++n ) {
            for ( idx_t v = 0; v < array.shape( 1 ); ++v ) {
                array( begin + n, v ) = tmp( order[n], v );
            }
        }
        return field.name();
    }
};

// ------------------------------------------------------------------

template <typename Value>
std::string reorder_field_T( Field& field, const std::vector<idx_t>& order, idx_t begin, idx_t end ) {
    if ( field.rank() == 1 ) {
        return ReorderField<Value, 1>::apply( field, order, begin, end );
    }
    else if ( field.rank() == 2 ) {
        return ReorderField<Value, 2>::apply( field, order, begin, end );
    }
    else {
        throw_Exception( "rank not supported", Here() );
    }
}

// ------------------------------------------------------------------

std::string reorder_field( Field& field, const std::vector<idx_t>& order, idx_t begin = 0,
                           idx_t end = std::numeric_limits<idx_t>::max() ) {
    if ( field.datatype() == array::DataType::kind<int>() ) {
        return reorder_field_T<int>( field, order, begin, end );
    }
    else if ( field.datatype() == array::DataType::kind<long>() ) {
        return reorder_field_T<long>( field, order, begin, end );
    }
    else if ( field.datatype() == array::DataType::kind<float>() ) {
        return reorder_field_T<float>( field, order, begin, end );
    }
    else if ( field.datatype() == array::DataType::kind<double>() ) {
        return reorder_field_T<double>( field, order, begin, end );
    }
    else {
        throw_Exception( "datatype not supported", Here() );
    }
}

// ------------------------------------------------------------------

void update_connectivity( mesh::HybridElements::Connectivity& connectivity, const std::vector<idx_t>& order ) {
    for ( idx_t b = 0; b < connectivity.blocks(); ++b ) {
        auto& block = connectivity.block( b );
        for ( idx_t r = 0; r < block.rows(); ++r ) {
            for ( idx_t c = 0; c < block.cols(); ++c ) {
                idx_t n = block( r, c );
                block.set( r, c, order.at( n ) );
            }
        }
    }
}

// ------------------------------------------------------------------

void reorder_connectivity( BlockConnectivity& connectivity, const std::vector<idx_t>& order ) {
    ATLAS_ASSERT( connectivity.rows() == static_cast<idx_t>( order.size() ) );
    BlockConnectivity tmp;
    tmp.add( connectivity.rows(), connectivity.cols(), connectivity.data(), true );
    for ( idx_t r = 0; r < connectivity.rows(); ++r ) {
        for ( idx_t c = 0; c < connectivity.cols(); ++c ) {
            if ( connectivity( r, c ) != tmp( r, c ) ) {
                ATLAS_DEBUG_VAR( r );
                ATLAS_DEBUG_VAR( c );
                ATLAS_DEBUG_VAR( connectivity( r, c ) );
                ATLAS_DEBUG_VAR( tmp( r, c ) );
            }
            ATLAS_ASSERT( connectivity( r, c ) == tmp( r, c ) );
            connectivity.set( r, c, tmp( order.at( r ), c ) );
        }
    }
}

// ------------------------------------------------------------------

void ReorderImpl::reorderNodes( Mesh& mesh, const std::vector<idx_t>& order ) {
    std::vector<idx_t> order_inverse( order.size() );
    for ( idx_t i = 0; i < static_cast<idx_t>( order.size() ); ++i ) {
        order_inverse[order[i]] = i;
    }

    for ( idx_t ifield = 0; ifield < mesh.nodes().nb_fields(); ++ifield ) {
        reorder_field( mesh.nodes().field( ifield ), order );
    }

    if ( mesh.cells().size() ) {
        update_connectivity( mesh.cells().node_connectivity(), order_inverse );
    }
    if ( mesh.edges().size() ) {
        update_connectivity( mesh.edges().node_connectivity(), order_inverse );
    }
}

// ------------------------------------------------------------------

void reorder_elements_using_nodes( Mesh& mesh, Mesh::HybridElements& elements ) {
    for ( idx_t t = 0; t < elements.nb_types(); ++t ) {
        auto& elems        = elements.elements( t );
        auto& connectivity = elems.node_connectivity();
        idx_t nb_nodes     = elems.nb_nodes();
        idx_t nb_elems     = elems.size();
        std::vector<std::pair<idx_t, idx_t>> node_lowest_index;
        node_lowest_index.reserve( elems.size() );
        for ( idx_t e = 0; e < nb_elems; ++e ) {
            idx_t lowest = std::numeric_limits<idx_t>::max();
            for ( idx_t n = 0; n < nb_nodes; ++n ) {
                lowest = std::min( lowest, connectivity( e, n ) );
            }
            node_lowest_index.emplace_back( lowest, e );
        }
        std::sort( node_lowest_index.begin(), node_lowest_index.end() );
        std::vector<idx_t> order;
        order.reserve( nb_elems );
        for ( const auto& pair : node_lowest_index ) {
            order.emplace_back( pair.second );
        }
        for ( idx_t ifield = 0; ifield < mesh.edges().nb_fields(); ++ifield ) {
            reorder_field( elements.field( ifield ), order, elems.begin(), elems.end() );
        }

        reorder_connectivity( elems.node_connectivity(), order );
    }
}

// ------------------------------------------------------------------

void ReorderImpl::reorderCellsUsingNodes( Mesh& mesh ) {
    reorder_elements_using_nodes( mesh, mesh.cells() );
}

// ------------------------------------------------------------------

void ReorderImpl::reorderEdgesUsingNodes( Mesh& mesh ) {
    reorder_elements_using_nodes( mesh, mesh.edges() );
}

// ------------------------------------------------------------------

void ReorderImpl::operator()( Mesh& mesh ) {
    ATLAS_TRACE( "ReorderImpl(mesh)" );

    reorderNodes( mesh, computeNodesOrder( mesh ) );
    reorderCellsUsingNodes( mesh );
    reorderEdgesUsingNodes( mesh );
}

// ------------------------------------------------------------------

namespace {
static ReorderBuilder<NoReorder> __ReorderReverseCuthillMckee( "none" );
}  // namespace

// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
