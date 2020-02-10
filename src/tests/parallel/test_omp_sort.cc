/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>   // generate, is_sorted
#include <functional>  // bind
#include <random>      // mt19937 and uniform_int_distribution
#include <vector>      // vector

#include "atlas/parallel/omp/sort.h"
#include "atlas/util/vector.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {


atlas::vector<int> create_random_data( int n ) {
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 eng( seed );  // a source of random data

    std::uniform_int_distribution<int> dist;
    atlas::vector<int> v( n );

    std::generate( v.begin(), v.end(), std::bind( dist, eng ) );
    return v;
}

//-----------------------------------------------------------------------------

CASE( "test_sort_little" ) {
    // Because there is little work to be done, the underlying algorithm will not use OpenMP here.
    auto integers = create_random_data( 20 );

    SECTION( "default order" ) {
        omp::sort( integers.begin(), integers.end() );
        EXPECT( std::is_sorted( integers.begin(), integers.end() ) );
    }

    SECTION( "reverse order" ) {
        omp::sort( integers.begin(), integers.end(), []( const int& a, const int& b ) { return b > a; } );
        EXPECT(
            std::is_sorted( integers.begin(), integers.end(), []( const int& a, const int& b ) { return b > a; } ) );
    }
}

CASE( "test_sort_large" ) {
    // There is enough work to use OpenMP underneith, if ATLAS_HAVE_OMP && OMP_NUM_THREADS > 1
    auto integers = create_random_data( 1000000 );

    SECTION( "default order" ) {
        omp::sort( integers.begin(), integers.end() );
        EXPECT( std::is_sorted( integers.begin(), integers.end() ) );
    }

    SECTION( "reverse order" ) {
        omp::sort( integers.begin(), integers.end(), []( const int& a, const int& b ) { return b > a; } );
        EXPECT(
            std::is_sorted( integers.begin(), integers.end(), []( const int& a, const int& b ) { return b > a; } ) );
    }
}

CASE( "test_merge_blocks" ) {
    // There is enough work to use OpenMP underneith, if ATLAS_HAVE_OMP && OMP_NUM_THREADS > 1

    auto nb_blocks  = 16;
    auto block_size = 100000;
    auto integers   = create_random_data( nb_blocks * block_size );
    std::vector<int> blocks_size( nb_blocks, block_size );

    auto end = integers.begin();
    for ( int b = 0; b < nb_blocks; ++b ) {
        auto begin = end;
        end        = begin + blocks_size[b];
        omp::sort( begin, end );
    }

    omp::merge_blocks( integers.begin(), integers.end(), blocks_size.begin(), blocks_size.end() );
    EXPECT( std::is_sorted( integers.begin(), integers.end() ) );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
