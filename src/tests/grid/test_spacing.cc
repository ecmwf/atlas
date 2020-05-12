/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <iomanip>
#include <sstream>

#include "eckit/types/Fraction.h"

#include "atlas/grid/Spacing.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/runtime/Log.h"

#include "tests/AtlasTestEnvironment.h"

using LinearSpacing = atlas::grid::LinearSpacing;
using Spacing       = atlas::grid::Spacing;
using Config        = atlas::util::Config;

namespace atlas {
namespace test {

using std::to_string;
std::string to_string( const eckit::Fraction& f ) {
    std::stringstream s;
    s << f;
    return s.str();
}

//-----------------------------------------------------------------------------

CASE( "LinearSpacing " ) {
    auto check = []( const Spacing& s, const std::vector<double>& expected ) {
        size_t j = 0;
        for ( double x : s ) {
            if ( not is_approximately_equal( x, expected[j++] ) ) {
                return false;
            }
        }
        return true;
    };

    /// Using the constructor LinearSpacing( start, end, N, endpoint ) we can create
    EXPECT( check( LinearSpacing( 2, 3, 5, true ), {2.0, 2.25, 2.5, 2.75, 3.0} ) );
    EXPECT( check( LinearSpacing( 2, 3, 5, false ), {2.0, 2.2, 2.4, 2.6, 2.8} ) );

    /// Using the constructor LinearSpacing( {start, end}, N, endpoint ) we can create
    EXPECT( check( LinearSpacing( {2, 3}, 5, true ), {2.0, 2.25, 2.5, 2.75, 3.0} ) );
    EXPECT( check( LinearSpacing( {2, 3}, 5, false ), {2.0, 2.2, 2.4, 2.6, 2.8} ) );

    /// Configuration parameters can be passed as well with following keys:
    EXPECT( check( Spacing( Config( "type", "linear" )( "start", 2 )( "end", 3 )( "N", 5 )( "endpoint", true ) ),
                   {2.0, 2.25, 2.5, 2.75, 3.0} ) );
    EXPECT( check( Spacing( Config( "type", "linear" )( "start", 2 )( "end", 3 )( "N", 5 )( "endpoint", false ) ),
                   {2.0, 2.2, 2.4, 2.6, 2.8} ) );

    /// Instead of the "end" key, you can provide the "length" key, to achieve the same results:
    EXPECT( check( Spacing( Config( "type", "linear" )( "start", 2 )( "length", 1 )( "N", 5 )( "endpoint", true ) ),
                   {2.0, 2.25, 2.5, 2.75, 3.0} ) );
    EXPECT( check( Spacing( Config( "type", "linear" )( "start", 2 )( "length", 1 )( "N", 5 )( "endpoint", false ) ),
                   {2.0, 2.2, 2.4, 2.6, 2.8} ) );
}

CASE( "LinearSpacing exactness ATLAS-276" ) {
    using eckit::Fraction;

    double north = 90.;
    double south = -90.;
    long N       = 217;

    SECTION( "LinearSpacing" ) {
        auto y = grid::LinearSpacing( north, south, N );

        // Ensure y does not go outside interval [90.,-90]
        EXPECT( !( y.front() > 90. ) );
        EXPECT( !( y.back() < -90. ) );

        // Exact front and end
        EXPECT( y.front() == 90. );
        EXPECT( y.back() == -90. );
    }

    SECTION( "LinearSpacing constructed with Fractions" ) {
        Fraction n( north );
        Fraction s( south );
        Fraction sn( 5, 6 );

        Fraction N_as_fraction = ( ( n - s ) / Fraction{5, 6} ) + 1;
        EXPECT( N_as_fraction.integer() );
        EXPECT_EQ( N_as_fraction.integralPart(), N );
        auto y = grid::LinearSpacing( n, s, N );

        EXPECT_EQ( Fraction{y.step()}, ( -Fraction{5, 6} ) );

        // Ensure y does not go outside interval [90.,-90]
        EXPECT( !( y.front() > 90. ) );
        EXPECT( !( y.back() < -90. ) );

        // Exact front and end
        EXPECT( y.front() == 90. );
        EXPECT( y.back() == -90. );
    }

    SECTION( "YSpace " ) {
        StructuredGrid::YSpace y( grid::LinearSpacing( north, south, N ) );

        // Ensure yspace does not go outside interval [90.,-90]
        EXPECT( !( y.front() > 90. ) );
        EXPECT( !( y.back() < -90. ) );

        // Exact front and end
        EXPECT( y.front() == 90. );
        EXPECT( y.back() == -90. );
    }

    SECTION( "LinearSpacing without endpoint" ) {
        auto y = grid::LinearSpacing( north, south, N - 1, /*endpoint*/ false );
        EXPECT_EQ( y.front(), 90. );
        EXPECT_EQ( Fraction{y.front() + ( N - 1 ) * y.step()}, -90. );
        EXPECT_EQ( Fraction{y.back() + y.step()}, -90. );
        EXPECT_EQ( Fraction{y.step()}, ( -Fraction{5, 6} ) );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
