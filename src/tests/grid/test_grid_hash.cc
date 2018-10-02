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

#include "atlas/domain.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"
#include "eckit/utils/MD5.h"

using Hash = eckit::MD5;

namespace atlas {
namespace test {

//----------------------------------------------------------------------------------------------------------------------

CASE( "test_global" ) {

    auto grids = {
      std::make_tuple( Grid( "O32" ),                                          "79599ee4825339e0ebbcd427d79d0177" ),
      std::make_tuple( Grid( "O32", Domain() ),                                "79599ee4825339e0ebbcd427d79d0177" ),
      std::make_tuple( Grid( "O32", RectangularDomain( {0,360}, {90,-90} ) ),  "79599ee4825339e0ebbcd427d79d0177" ),
      std::make_tuple( Grid( "O32", ZonalBandDomain( {90,-90} ) ),             "79599ee4825339e0ebbcd427d79d0177" ),
    };

    int c{0};
    for( auto entry: grids ) {
        Grid grid = std::get<0>( entry );
        std::string hash = std::get<1>( entry );
        SECTION( "grid: " + std::to_string(c) ) {
            Hash h; grid.hash( h );
            if( hash.empty() ) {
                Log::info() << "grid " << c << "   hash = " << std::string(h) << std::endl;
            }
            EXPECT( std::string(h) == hash );
        }
        c++;
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
