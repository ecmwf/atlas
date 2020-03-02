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

#include "atlas/grid.h"
#include "atlas/parallel/omp/copy.h"
#include "atlas/util/Point.h"
#include "atlas/util/vector.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test omp::copy" ) {
    Grid grid( "O400" );
    atlas::vector<PointXY> points( grid.size() );
    omp::copy( grid.xy().begin(), grid.xy().end(), points.begin() );
}


CASE( "test atlas::vector assign" ) {
    // atlas::vector assign uses omp::copy
    // It should reproduce exactly the same points as previous test case
    Grid grid( "O400" );
    atlas::vector<PointXY> points;
    points.assign( grid.xy().begin(), grid.xy().end() );
    EXPECT( points.size() == grid.size() );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
