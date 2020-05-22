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

#include "atlas/library/config.h"
#include "atlas/util/vector.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------
//


template <typename T>
atlas::vector<T> square( const int n ) {
    atlas::vector<T> x( n );
    for ( int i = 0; i < n; i++ )
        x[i] = static_cast<T>( i ) * static_cast<T>( i );

    return x;
}

template <typename T>
void pp( const std::string& t, T& v ) {
    Log::info() << t << " = " << std::endl;
    Log::info() << v.size() << std::endl;
    for ( int i = 0; i < v.size(); i++ )
        Log::info() << i << " > " << v[i] << std::endl;
}


CASE( "test_vector" ) {
    const int N = 100;

    // Ctor
    atlas::vector<double> x( N );

    std::iota( std::begin( x ), std::end( x ), 0 );

    // Copy ctor
    atlas::vector<double> y = x;

    EXPECT( x.size() == y.size() );

    for ( int i = 0; i < x.size(); i++ ) {
        EXPECT( x[i] == y[i] );
    }


    // Assignment operator
    auto z = square<int>( 20 );

    for ( int i = 0; i < z.size(); i++ ) {
        EXPECT( z[i] == i * i );
    }

    // Assignment operator
    x = y;

    for ( int i = 0; i < x.size(); i++ ) {
        EXPECT( static_cast<int>( x[i] ) == i );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
