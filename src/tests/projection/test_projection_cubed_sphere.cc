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
    Grid gEA{ "CS-EA-L-" + std::to_string(resolution)};
    Grid gLFR{ "CS-LFR-L-" + std::to_string(resolution)};

    using util::Constants;

    util::Config params;
    CubedSphereTiles f("cubedsphere_fv3");
    CubedSphereTiles l("cubedsphere_lfric");

    double cd[2];

    idx_t jn(0);

    std::array<idx_t, 7> EAOffset{0,
                                  resolution * resolution + 1,
                                  2 * resolution * resolution + 2,
                                  3 * resolution * resolution + 2,
                                  4 * resolution * resolution + 2,
                                  5 * resolution * resolution + 2,
                                  6 * resolution * resolution + 2};

    for ( auto crd : gEA.lonlat() ) {
       atlas::PointLonLat pointLonLat = crd;
       cd[LON] = pointLonLat.lon();
       cd[LAT] = pointLonLat.lat();

       int t = f.tileFromLonLat(cd);  

       gEA.projection().lonlat2xy(crd);
       cd[LON] = crd.lon();
       cd[LAT] = crd.lat();

       int t2 = f.tileFromXY(cd);

       for (std::size_t i = 0; i < 6; ++i) {
           if (jn >= EAOffset[i] && jn < EAOffset[i+1]) {
               EXPECT(t == static_cast<idx_t>(i));
               EXPECT(t2 == static_cast<idx_t>(i));
           }
       }
       ++jn;
    }

    std::array<idx_t, 7> LFRicOffset{0,
                                  resolution * resolution,
                                  2 * resolution * resolution,
                                  3 * resolution * resolution,
                                  4 * resolution * resolution,
                                  4 * resolution * resolution + (resolution + 1) * (resolution + 1),
                                  6 * resolution * resolution + 2};

    jn = 0;
    for ( auto crd : gLFR.lonlat() ) {
       atlas::PointLonLat pointLonLat = crd;
       cd[LON] = pointLonLat.lon();
       cd[LAT] = pointLonLat.lat();

       int t3 = l.tileFromLonLat(cd);

       gLFR.projection().lonlat2xy(crd);
       cd[LON] = crd.lon();
       cd[LAT] = crd.lat();

       int t4 = l.tileFromXY(cd);

       for (std::size_t i = 0; i < 6; ++i) {
           if (jn >= LFRicOffset[i] && jn < LFRicOffset[i+1]) {
               EXPECT(t3 == static_cast<idx_t>(i));
               EXPECT(t4 == static_cast<idx_t>(i));
           }
       }
       ++jn;
    }

}

CASE( "test_projection_cubedsphere_xy_latlon" ) {
    int resolution( 12 );
    std::vector<std::string> grid_names{"CS-EA-L-" + std::to_string( resolution ),
                                        "CS-ED-L-" + std::to_string( resolution )};

    for ( std::string& s : grid_names ) {
        Grid g{s};
        for ( auto crd : g.lonlat() ) {
            Point2 lonlat{crd};
            g->projection().lonlat2xy( crd );
            g->projection().xy2lonlat( crd );
            // except for point lonlat (90,82.5) on compiler pgc++
            // we have a maximum error tolerance of 1e-11
            EXPECT_APPROX_EQ( lonlat, crd, 1e-6 );
        }
        for ( auto crd : g.xy() ) {
            Point2 xy{crd};
            g->projection().xy2lonlat( crd );
            g->projection().lonlat2xy( crd );
            EXPECT_APPROX_EQ( xy, crd, 1e-6 );
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
