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
#include <limits>
#include <vector>

#include "atlas/grid.h"
#include "atlas/util/KDTree.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::util::KDTree;

namespace atlas {
namespace test {

CASE( "test kdtree" ) {
    auto grid = Grid{"O32"};

    KDTree<idx_t> search;
    search.reserve( grid.size() );
    idx_t n{0};
    for ( auto& point : grid.lonlat() ) {
        search.insert( point, n++ );
    }
    search.build();
    EXPECT_NO_THROW( search.nearestNeighbour( PointLonLat{180., 45.} ) );
    // ...
    // Search 4 nearest neighbours (k=4), sorted by shortest distance
    auto neighbours          = search.kNearestNeighbours( PointLonLat{180., 45.}, 4 ).payloads();
    auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761};
    EXPECT( neighbours == expected_neighbours );
}

CASE( "test assertion" ) {
    auto grid = Grid{"O32"};

    KDTree<idx_t> search;
    search.reserve( grid.size() );
    idx_t n{0};
    for ( auto& point : grid.lonlat() ) {
        search.insert( point, n++ );
    }
    // Forgot to call search.build() --> assertion thrown when trying to access
    EXPECT_THROWS_AS( search.nearestNeighbour( PointLonLat{180., 45.} ), eckit::AssertionFailed );
}

CASE( "test no assertion" ) {
    // Like case "test assertion", but without reserving size
    auto grid = Grid{"O32"};

    KDTree<idx_t> search;
    // No search.reserve() --> build() will not be necessary.
    idx_t n{0};
    for ( auto& point : grid.lonlat() ) {
        search.insert( point, n++ );
    }
    // search.build() Not required
    EXPECT_NO_THROW( search.nearestNeighbour( PointLonLat{180., 45.} ) );
}

const KDTree<idx_t>& search() {
    static KDTree<idx_t>* kdtree = []() {
        KDTree<idx_t>* kdtree = new KDTree<idx_t>;
        auto grid             = Grid{"O32"};
        kdtree->reserve( grid.size() );
        idx_t n{0};
        for ( auto& point : grid.lonlat() ) {
            kdtree->insert( point, n++ );
        }
        kdtree->build();
        return kdtree;
    }();
    EXPECT( kdtree );
    return *kdtree;
}

CASE( "test nearestNeighbour" ) {
    auto neighbour          = search().nearestNeighbour( PointLonLat{180., 45.} ).payload();
    auto expected_neighbour = 760;
    EXPECT( neighbour == expected_neighbour );
}

CASE( "test kNearestNeighbours" ) {
    auto neighbours          = search().kNearestNeighbours( PointLonLat{180., 45.}, 4 ).payloads();
    auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761};
    EXPECT( neighbours == expected_neighbours );
}

CASE( "test findInSphere" ) {
    constexpr double km      = 1000.;
    auto neighbours          = search().findInSphere( PointLonLat{180., 45.}, 500 * km ).payloads();
    auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761, 841, 843, 682};
    EXPECT( neighbours == expected_neighbours );
    Log::info() << neighbours << std::endl;
}

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
