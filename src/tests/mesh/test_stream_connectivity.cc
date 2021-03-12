/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/io/Buffer.h"
#include "eckit/serialisation/ResizableMemoryStream.h"

#include "atlas/mesh/Connectivity.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::mesh;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_stream_irregular_connectivity" ) {
    eckit::Buffer b{0};
    eckit::ResizableMemoryStream s{b};

    // Create stream
    {
        IrregularConnectivity conn( "bla" );

        constexpr idx_t vals[4] = {2, 3, 5, 6};
        conn.add( 1, 4, vals, /* fortran-array = */ false );

        constexpr idx_t vals2[6] = {1, 3, 4, 3, 7, 8};
        conn.add( 2, 3, vals2, /* fortran-array = */ false );
        s << conn;
    }


    s.rewind();
    // Read from stream
    {
        IrregularConnectivity conn( s );
        EXPECT( conn.rows() == 3 );
        EXPECT( conn.cols( 0 ) == 4 );
        EXPECT( conn.cols( 1 ) == 3 );
        EXPECT( conn.cols( 2 ) == 3 );
        EXPECT( conn( 0, 0 ) == 2 );
        EXPECT( conn( 0, 1 ) == 3 );
        EXPECT( conn( 0, 2 ) == 5 );
        EXPECT( conn( 0, 3 ) == 6 );
        EXPECT( conn( 1, 0 ) == 1 );
        EXPECT( conn( 1, 1 ) == 3 );
        EXPECT( conn( 1, 2 ) == 4 );
        EXPECT( conn( 2, 0 ) == 3 );
        EXPECT( conn( 2, 1 ) == 7 );
        EXPECT( conn( 2, 2 ) == 8 );
        EXPECT( conn.name() == "bla" );
    }


    s.rewind();
    // Read from stream
    {
        IrregularConnectivity conn;
        s >> conn;
        EXPECT( conn.rows() == 3 );
        EXPECT( conn.cols( 0 ) == 4 );
        EXPECT( conn.cols( 1 ) == 3 );
        EXPECT( conn.cols( 2 ) == 3 );
        EXPECT( conn( 0, 0 ) == 2 );
        EXPECT( conn( 0, 1 ) == 3 );
        EXPECT( conn( 0, 2 ) == 5 );
        EXPECT( conn( 0, 3 ) == 6 );
        EXPECT( conn( 1, 0 ) == 1 );
        EXPECT( conn( 1, 1 ) == 3 );
        EXPECT( conn( 1, 2 ) == 4 );
        EXPECT( conn( 2, 0 ) == 3 );
        EXPECT( conn( 2, 1 ) == 7 );
        EXPECT( conn( 2, 2 ) == 8 );
        EXPECT( conn.name() == "bla" );
    }
}

//-----------------------------------------------------------------------------

CASE( "test_stream_block_connectivity" ) {
    eckit::Buffer b{0};
    eckit::ResizableMemoryStream s{b};

    // Create stream
    {
        BlockConnectivity conn( 2, 3, {1, 3, 4, 3, 7, 8} );
        s << conn;
    }


    s.rewind();
    // Read from stream
    {
        BlockConnectivity conn( s );
        EXPECT( conn.rows() == 2 );
        EXPECT( conn.cols() == 3 );
        EXPECT( conn( 0, 0 ) == 1 );
        EXPECT( conn( 0, 1 ) == 3 );
        EXPECT( conn( 0, 2 ) == 4 );
        EXPECT( conn( 1, 0 ) == 3 );
        EXPECT( conn( 1, 1 ) == 7 );
        EXPECT( conn( 1, 2 ) == 8 );
    }


    s.rewind();
    // Read from stream
    {
        BlockConnectivity conn;
        s >> conn;
        EXPECT( conn.rows() == 2 );
        EXPECT( conn.cols() == 3 );
        EXPECT( conn( 0, 0 ) == 1 );
        EXPECT( conn( 0, 1 ) == 3 );
        EXPECT( conn( 0, 2 ) == 4 );
        EXPECT( conn( 1, 0 ) == 3 );
        EXPECT( conn( 1, 1 ) == 7 );
        EXPECT( conn( 1, 2 ) == 8 );
    }
}

CASE( "test_stream_multiblock_connectivity" ) {
    eckit::Buffer b{0};
    eckit::ResizableMemoryStream s{b};

    // Create stream
    {
        MultiBlockConnectivity conn( "mbc" );
        conn.add( BlockConnectivity{3, 5, {3, 7, 1, 4, 5, 6, 4, 56, 8, 4, 1, 3, 76, 4, 3}} );
        conn.add( BlockConnectivity{3, 2, {4, 75, 65, 45, 51, 35}} );
        s << conn;
    }

    s.rewind();
    // Read from stream
    {
        MultiBlockConnectivity conn{s};
        EXPECT( conn.name() == "mbc" );
        EXPECT( conn.blocks() == 2 );
        EXPECT( conn( 3, 1 ) == 75 );
        EXPECT( conn( 4, 1 ) == 45 );
        EXPECT( conn( 5, 0 ) == 51 );
        EXPECT( conn( 1, 0, 1 ) == 75 );
        EXPECT( conn( 1, 1, 1 ) == 45 );
        EXPECT( conn( 1, 2, 0 ) == 51 );
    }

    s.rewind();
    // Read from stream
    {
        MultiBlockConnectivity conn;
        s >> conn;
        EXPECT( conn.name() == "mbc" );
        EXPECT( conn.blocks() == 2 );
        EXPECT( conn( 3, 1 ) == 75 );
        EXPECT( conn( 4, 1 ) == 45 );
        EXPECT( conn( 5, 0 ) == 51 );
        EXPECT( conn( 1, 0, 1 ) == 75 );
        EXPECT( conn( 1, 1, 1 ) == 45 );
        EXPECT( conn( 1, 2, 0 ) == 51 );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
