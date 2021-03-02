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
#include <list>
#include <queue>
#include <utility>
#include <vector>


#include "atlas/mesh/actions/ReorderReverseCuthillMckee.h"

#include "atlas/array.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/detail/AccumulateFacets.h"
#include "atlas/runtime/Trace.h"

#include "eckit/linalg/SparseMatrix.h"
#include "eckit/linalg/Triplet.h"

namespace atlas {
namespace mesh {
namespace actions {

using NotVisited_t = std::list<std::pair<int, double>>;

NotVisited_t::iterator find_index( NotVisited_t& a, int x ) {
    return std::find_if( a.begin(), a.end(), [x]( NotVisited_t::value_type& element ) { return element.first == x; } );
}

class CuthillMckee {
private:
    eckit::linalg::SparseMatrix sparse_;

public:
    CuthillMckee( eckit::linalg::SparseMatrix&& m ) : sparse_( std::move( m ) ) {}

    // Function to generate degree of all the nodes
    std::vector<int> computeDegrees( const eckit::linalg::SparseMatrix& m ) {
        std::vector<int> degrees;
        degrees.reserve( m.rows() );
        const auto outer = m.outer();
        for ( size_t i = 0; i < m.rows(); i++ ) {
            degrees.emplace_back( outer[i + 1] - outer[i] );
        }
        return degrees;
    }

    // Cuthill-Mckee algorithm
    std::vector<idx_t> order() {
        auto degrees = computeDegrees( sparse_ );

        std::queue<idx_t> Q;
        std::vector<idx_t> R;
        R.reserve( degrees.size() );

        std::vector<idx_t> sorted_row;
        sorted_row.reserve( sparse_.cols() );

        NotVisited_t not_visited;

        for ( size_t i = 0; i < degrees.size(); i++ ) {
            not_visited.emplace_back( i, degrees[i] );
        }

        not_visited.sort( []( const NotVisited_t::value_type& a, const NotVisited_t::value_type& b ) {
            return a.second < b.second;
        } );

        // notVisited helps in running Breadth First Search even when there are disjoint graphs

        const auto outer = sparse_.outer();
        const auto index = sparse_.inner();

        while ( not_visited.size() ) {
            Q.emplace( not_visited.front().first );
            not_visited.pop_front();

            // Simple Breadth First Search
            while ( !Q.empty() ) {
                idx_t row = Q.front();

                sorted_row.clear();

                for ( auto j = outer[row]; j < outer[row + 1]; ++j ) {
                    if ( index[j] != row ) {
                        auto found_index = find_index( not_visited, index[j] );
                        if ( found_index != not_visited.end() ) {
                            sorted_row.emplace_back( index[j] );
                            not_visited.erase( found_index );
                        }
                    }
                }

                std::sort( sorted_row.begin(), sorted_row.end(),
                           [&]( idx_t i, idx_t j ) { return degrees[i] - degrees[j]; } );

                for ( size_t i = 0; i < sorted_row.size(); i++ ) {
                    Q.emplace( sorted_row[i] );
                }

                R.emplace_back( Q.front() );
                Q.pop();
            }
        }

        return R;
    }
};

class ReverseCuthillMckee {
private:
    CuthillMckee cuthill_mckee_;

public:
    ReverseCuthillMckee( eckit::linalg::SparseMatrix&& m ) : cuthill_mckee_( std::move( m ) ) {}

    // Reverse Cuthill-Mckee algorithm
    std::vector<idx_t> order() {
        std::vector<idx_t> cuthill = cuthill_mckee_.order();

        idx_t size = static_cast<idx_t>( cuthill.size() );
        idx_t n    = size;

        if ( n % 2 == 0 ) {
            n -= 1;
        }

        n = n / 2;

        for ( idx_t i = 0; i <= n; i++ ) {
            idx_t j               = cuthill[size - 1 - i];
            cuthill[size - 1 - i] = cuthill[i];
            cuthill[i]            = j;
        }

        return cuthill;
    }
};


ReorderReverseCuthillMckee::ReorderReverseCuthillMckee( const eckit::Parametrisation& config ) {
    config.get( "ghost_at_end", ghost_at_end_ );
}

std::vector<idx_t> ReorderReverseCuthillMckee::computeNodesOrder( Mesh& mesh ) {
    ATLAS_TRACE();
    std::vector<eckit::linalg::Triplet> n2n_triplets;
    {
        auto ghost_node = array::make_view<int, 1>( mesh.nodes().ghost() );

        std::vector<idx_t> facet_nodes_data;
        std::vector<idx_t> connectivity_facet_to_elem;
        idx_t nb_facets;
        idx_t nb_inner_facets;
        idx_t missing_value;

        detail::accumulate_facets( mesh.cells(), mesh.nodes(),
                                   facet_nodes_data,  // shape(nb_facets,nb_nodes_per_facet)
                                   connectivity_facet_to_elem, nb_facets, nb_inner_facets, missing_value );

        n2n_triplets.reserve( 2 * nb_facets );

        for ( idx_t e = 0; e < nb_facets; ++e ) {
            idx_t node_0 = facet_nodes_data[2 * e + 0];
            idx_t node_1 = facet_nodes_data[2 * e + 1];
            if ( ghost_at_end_ && ( ghost_node( node_0 ) || ghost_node( node_1 ) ) ) {
                continue;
            }
            n2n_triplets.emplace_back( node_0, node_1, 1. );
            n2n_triplets.emplace_back( node_1, node_0, 1. );
        }
    }

    std::sort( n2n_triplets.begin(), n2n_triplets.end() );

    auto nb_nodes = static_cast<eckit::linalg::Size>( mesh.nodes().size() );
    return ReverseCuthillMckee{{nb_nodes, nb_nodes, n2n_triplets}}.order();
}


namespace {
static ReorderBuilder<ReorderReverseCuthillMckee> __ReorderReverseCuthillMckee( "reverse_cuthill_mckee" );
}  // namespace

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
