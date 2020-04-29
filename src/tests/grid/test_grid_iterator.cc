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

#include "eckit/types/FloatCompare.h"
#include "eckit/types/Fraction.h"

#include "atlas/grid/Iterator.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

bool operator==( const std::vector<PointLonLat>& ll, const std::vector<PointXY>& xy ) {
    if ( ll.size() != xy.size() ) {
        return false;
    }
    for ( size_t i = 0; i < xy.size(); ++i ) {
        if ( ll[i] != xy[i] ) {
            return false;
        }
    }
    return true;
}

//-----------------------------------------------------------------------------

CASE( "test_iterator" ) {
    std::vector<Grid> grids;

    std::vector<PointXY> points{{0, 90},  {90, 90}, {180, 90}, {270, 90}, {0, 0},     {90, 0},
                                {180, 0}, {270, 0}, {0, -90},  {90, -90}, {180, -90}, {270, -90}};

    grids.emplace_back( "L4x3" );
    grids.emplace_back( UnstructuredGrid{points} );

    for ( auto grid : grids ) {
        Log::debug() << "grid : " << grid.name() << std::endl;

        std::vector<PointXY> points_xy;
        std::vector<PointLonLat> points_lonlat;

        // Iteration of xy in range-based for
        for ( const PointXY& xy : grid.xy() ) {
            points_xy.push_back( xy );
        }
        EXPECT( points_xy.size() == static_cast<size_t>( grid.size() ) );
        EXPECT( points_xy == points );

        // Iteration of lonlat in range-based for
        for ( const PointLonLat& ll : grid.lonlat() ) {
            points_lonlat.push_back( ll );
        }
        EXPECT( points_lonlat.size() == static_cast<size_t>( grid.size() ) );
        EXPECT( points_lonlat == points );


        // Use of xy iterator in STL algorithms
        std::copy( grid.xy().begin(), grid.xy().end(), points_xy.begin() );
        EXPECT( points_xy == points );

        // Use of lonlat iterator in STL algorithms
        std::copy( grid.lonlat().begin(), grid.lonlat().end(), points_lonlat.begin() );
        EXPECT( points_lonlat == points );

        // Random access for xy iterator, required for parallelisation
        EXPECT( *( grid.xy().begin() + grid.size() / 2 ) == PointXY( 180., 0. ) );
        EXPECT( grid.xy().begin() + grid.size() == grid.xy().end() );

        // Random access for lonlat iterator, required for parallelisation
        EXPECT( *( grid.lonlat().begin() + grid.size() / 2 ) == PointLonLat( 180., 0. ) );
        EXPECT( grid.lonlat().begin() + grid.size() == grid.lonlat().end() );
    }
}

//-----------------------------------------------------------------------------

CASE( "ATLAS-276" ) {
    using eckit::Fraction;
    using eckit::types::is_approximately_equal;

    Fraction n( 90. );
    Fraction s( -90. );
    Fraction sn( 5, 6 );

    Fraction N_as_fraction = ( ( n - s ) / sn ) + 1;
    EXPECT( N_as_fraction.integer() );

    long N = N_as_fraction.integralPart();
    EXPECT( N == 217 );

    StructuredGrid::YSpace y( grid::LinearSpacing( n, s, N ) );

    // Tolerance might be suitable
    EXPECT( is_approximately_equal( y.front(), 90. ) );
    EXPECT( is_approximately_equal( y.back(), -90. ) );

    // Tolerance isn't suitable
    EXPECT( !( y.front() > 90. ) );
    EXPECT( !( y.back() < -90. ) );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
