/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/option.h"
#include "atlas/array.h"
#include "atlas/grid/Vertical.h"
#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::util;


namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

std::vector<double> IFS_vertical_coordinates( idx_t nlev ) {
    std::vector<double> zcoord( nlev + 2 );
    zcoord[0]        = 0.;
    zcoord[nlev + 1] = 1.;
    double dzcoord   = 1. / double( nlev );
    for ( idx_t jlev = 1; jlev <= nlev; ++jlev ) {
        zcoord[jlev] = jlev * dzcoord - 0.5 * dzcoord;
    }
    return zcoord;
}

std::vector<double> zrange( idx_t nlev, double min, double max ) {
    std::vector<double> zcoord( nlev );
    double dzcoord = ( max - min ) / double( nlev - 1 );
    for ( idx_t jlev = 0; jlev < nlev; ++jlev ) {
        zcoord[jlev] = min + jlev * dzcoord;
    }
    return zcoord;
}


CASE( "test vertical;  default" ) {
    Vertical vertical;
    EXPECT( vertical.size() == 0 );
    EXPECT( vertical.k_begin() == 0 );
    EXPECT( vertical.k_end() == 0 );
    EXPECT( vertical.boundaries() == false );
}
CASE( "test vertical;  config levels without boundaries" ) {
    Vertical vertical( option::levels( 10 ) );
    EXPECT( vertical.size() == 10 );
    EXPECT( vertical.k_begin() == 0 );
    EXPECT( vertical.k_end() == 10 );
    EXPECT( vertical.boundaries() == false );
}
CASE( "test vertical;  config levels with boundaries" ) {
    Vertical vertical( option::levels( 10 ) | Config( "boundaries", true ) );
    EXPECT( vertical.size() == 12 );
    EXPECT( vertical.k_begin() == 1 );
    EXPECT( vertical.k_end() == 11 );
    EXPECT( vertical.boundaries() == true );
}
CASE( "test vertical;  array with boundaries" ) {
    Vertical vertical( 5, IFS_vertical_coordinates(5), Config("boundaries",true ) );
    EXPECT( vertical.size() == 7 );
    EXPECT( vertical.k_begin() == 1 );
    EXPECT( vertical.k_end() == 6 );
    EXPECT( vertical.boundaries() == true );
}
CASE( "test vertical;  array without boundaries" ) {
    Vertical vertical( 5, zrange(5, 0., 1.) );
    EXPECT( vertical.size() == 5 );
    EXPECT( vertical.k_begin() == 0 );
    EXPECT( vertical.k_end() == 5 );
    EXPECT( vertical.boundaries() == false );
}

//-----------------------------------------------------------------------------


}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
