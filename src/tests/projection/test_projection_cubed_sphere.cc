/*
 * (C) Crown Copyright 2021 MetOffice
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <vector>
#include "atlas/grid.h"
#include "atlas/grid/Grid.h"
#include "atlas/projection/Projection.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {


//-----------------------------------------------------------------------------

CASE( "test_projection_cubedsphere_xy_latlon" ) {

    int resolution(12);
    std::vector<std::string> grid_names{"CS-EA-" + std::to_string(resolution),
                                        "CS-ED-" + std::to_string(resolution)};

    for (std::string & s : grid_names) {
        Grid g{ s };
        for ( auto crd : g.lonlat() ) {
            Point2 lonlat { crd };
            g->projection().lonlat2xy( crd );
            g->projection().xy2lonlat( crd );
            // except for point lonlat (90,82.5) on compiler pgc++
            // we have a error of 1e-11
            EXPECT_APPROX_EQ(lonlat, crd, 1e-6);
        }
        for ( auto crd : g.xy() ) {
            Point2 xy { crd };
            g->projection().xy2lonlat( crd );
            g->projection().lonlat2xy( crd );
            EXPECT_APPROX_EQ(xy, crd, 1e-6);
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
