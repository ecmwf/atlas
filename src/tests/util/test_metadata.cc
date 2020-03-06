/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array/ArrayView.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Metadata.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_broadcast_to_self" ) {
    Metadata metadata;
    if ( mpi::comm().rank() == 0 ) {
        metadata.set( "paramID", 128 );
    }

    // broadcast
    metadata.broadcast();

    EXPECT( metadata.has( "paramID" ) );
    if ( metadata.has( "paramID" ) ) {
        EXPECT( metadata.get<int>( "paramID" ) == 128 );
    }
}

// -----------------------------------------------------------------------------

CASE( "test_broadcast_to_other" ) {
    size_t root = 0;
    Metadata global;
    if ( mpi::comm().rank() == root ) {
        global.set( "paramID", 128 );
    }

    Metadata local;

    // broadcast
    global.broadcast( local );

    EXPECT( local.has( "paramID" ) );
    if ( local.has( "paramID" ) ) {
        EXPECT( local.get<int>( "paramID" ) == 128 );
    }

    if ( mpi::comm().rank() != root ) {
        EXPECT( !global.has( "paramID" ) );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
