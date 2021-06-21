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
#include "atlas/util/Constants.h"


#include "atlas/grid/Tiles.h"
#include "atlas/grid/detail/tiles/Tiles.h"

#include "atlas/grid/detail/tiles/FV3Tiles.h"
#include "atlas/grid/detail/tiles/LFRicTiles.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {


//-----------------------------------------------------------------------------

CASE ("test_tiles") {
   int resolution(2);
    Grid g{ "CS-EA-" + std::to_string(resolution)};

    using atlas::cubedspheretiles::FV3CubedSphereTiles;
    using atlas::cubedspheretiles::LFRicCubedSphereTiles;
    using util::Constants;

    util::Config params;
    FV3CubedSphereTiles f(params);
    LFRicCubedSphereTiles l(params);

    double cd[2];
    for ( auto crd : g.lonlat() ) {
       std::cout << "fv3 crd "  << crd[LON] << " " << crd[LAT]  <<  std::endl;
       atlas::PointLonLat pointLonLat = crd * Constants::degreesToRadians();
       cd[LON] = pointLonLat.lon();
       cd[LAT] = pointLonLat.lat();
       std::cout << "fv3 cd(lon , lat) = "
                 << cd[LON] << " " << cd[LAT]  << std::endl;

       int t = f.tileFromLonLat(cd);
       int t2 = l.tileFromLonLat(cd);
       std::cout << "fv3 lfric tilesfromlonlat = " << t << " " << t2 << std::endl;

       g.projection().lonlat2xy(crd);
       cd[LON] = crd.lon();
       cd[LAT] = crd.lat();

       t = f.tileFromXY(cd);
       t2 = l.tileFromXY(cd);
       std::cout << "fv3 lfric tilesfromxy = " << t << " " << t2 << std::endl;
    }

    params.set("tile type", "LFRicCubedSphereTiles");
    atlas::CubedSphereTiles tiles(params);
    for ( auto crd : g.lonlat() ) {
        atlas::PointLonLat pointLonLat = crd * Constants::degreesToRadians();
        cd[LON] = pointLonLat.lon();
        cd[LAT] = pointLonLat.lat();
        std::cout << "generic cd(lon , lat) = "
                  << cd[LON] << " " << cd[LAT]  << std::endl;
        std::cout <<  f.tileFromLonLat(cd)   << std::endl;
        std::cout <<  tiles.tileFromLonLat(cd)   << std::endl;

    }



}

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
            // we have a maximum error tolerance of 1e-11
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
