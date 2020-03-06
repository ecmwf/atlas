/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/library/defines.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::mesh;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

#if ATLAS_HAVE_FORTRAN
#define FORTRAN_BASE 1
#define FROM_FORTRAN -1
#define TO_FORTRAN +1
#define IN_FORTRAN -1
#else
#define FORTRAN_BASE 0
#define FROM_FORTRAN
#define TO_FORTRAN
#define IN_FORTRAN
#endif

CASE( "test_irregular_connectivity" ) {
    IrregularConnectivity conn( "mesh" );
    EXPECT( conn.rows() == 0 );
    EXPECT( conn.maxcols() == 0 );

    constexpr idx_t vals[4] = {2, 3, 5, 6};
    conn.add( 1, 4, vals, /* fortran-array = */ false );

    EXPECT( conn.rows() == 1 );
    EXPECT( conn.cols( 0 ) == 4 );
    EXPECT( conn.mincols() == 4 );
    EXPECT( conn.maxcols() == 4 );

    EXPECT( conn( 0, 0 ) == 2 );
    EXPECT( conn( 0, 1 ) == 3 );
    EXPECT( conn( 0, 2 ) == 5 );
    EXPECT( conn( 0, 3 ) == 6 );

    EXPECT( conn.row( 0 )( 0 ) == 2 );
    EXPECT( conn.row( 0 )( 1 ) == 3 );
    EXPECT( conn.row( 0 )( 2 ) == 5 );
    EXPECT( conn.row( 0 )( 3 ) == 6 );

    constexpr idx_t vals2[6] = {1, 3, 4, 3, 7, 8};
    conn.add( 2, 3, vals2, /* fortran-array = */ true );

    conn.dump( Log::info() );
    Log::info() << std::endl;

    EXPECT( conn.rows() == 3 );
    EXPECT( conn.cols( 1 ) == 3 );
    EXPECT( conn.cols( 2 ) == 3 );
    EXPECT( conn.mincols() == 3 );
    EXPECT( conn.maxcols() == 4 );

    EXPECT( conn( 1, 0 ) == 1 IN_FORTRAN );
    EXPECT( conn( 1, 1 ) == 3 IN_FORTRAN );
    EXPECT( conn( 1, 2 ) == 4 IN_FORTRAN );

    EXPECT( conn.row( 2 )( 0 ) == 3 IN_FORTRAN );
    EXPECT( conn.row( 2 )( 1 ) == 7 IN_FORTRAN );
    EXPECT( conn.row( 2 )( 2 ) == 8 IN_FORTRAN );

    conn.set( 1, 1, 9 IN_FORTRAN );
    EXPECT( conn( 1, 0 ) == 1 IN_FORTRAN );
    EXPECT( conn( 1, 1 ) == 9 IN_FORTRAN );
    EXPECT( conn( 1, 2 ) == 4 IN_FORTRAN );

    constexpr idx_t vals3[3] = {6 IN_FORTRAN, 7 IN_FORTRAN, 5 IN_FORTRAN};
    conn.set( 2, vals3 );
    EXPECT( conn( 2, 0 ) == 6 IN_FORTRAN );
    EXPECT( conn( 2, 1 ) == 7 IN_FORTRAN );
    EXPECT( conn( 2, 2 ) == 5 IN_FORTRAN );

    constexpr idx_t vals4[8] = {2, 11, 51, 12, 4, 13, 55, 78};

    conn.insert( 1, 2, 4, vals4, /* fortran-array = */ false );
    EXPECT( conn.mincols() == 3 );
    EXPECT( conn.maxcols() == 4 );

    EXPECT( conn.rows() == 5 );
    EXPECT( conn.cols( 0 ) == 4 );
    EXPECT( conn.cols( 1 ) == 4 );
    EXPECT( conn.cols( 2 ) == 4 );
    EXPECT( conn.cols( 3 ) == 3 );
    EXPECT( conn.cols( 4 ) == 3 );

    EXPECT( conn( 1, 0 ) == 2 );
    EXPECT( conn( 1, 1 ) == 11 );
    EXPECT( conn( 1, 2 ) == 51 );
    EXPECT( conn( 1, 3 ) == 12 );

    EXPECT( conn( 2, 0 ) == 4 );
    EXPECT( conn( 2, 1 ) == 13 );
    EXPECT( conn( 2, 2 ) == 55 );
    EXPECT( conn( 2, 3 ) == 78 );

    EXPECT( conn( 3, 0 ) == 1 IN_FORTRAN );
    EXPECT( conn( 3, 1 ) == 9 IN_FORTRAN );
    EXPECT( conn( 3, 2 ) == 4 IN_FORTRAN );

    constexpr idx_t vals5[2] = {3, 6};
    conn.insert( 3, 1, 2, vals5, true );

    EXPECT( conn.mincols() == 2 );
    EXPECT( conn.maxcols() == 4 );

    EXPECT( conn.rows() == 6 );
    EXPECT( conn.cols( 3 ) == 2 );
    EXPECT( conn.cols( 4 ) == 3 );

    EXPECT( conn( 3, 0 ) == 3 IN_FORTRAN );
    EXPECT( conn( 3, 1 ) == 6 IN_FORTRAN );

    EXPECT( conn( 4, 0 ) == 1 IN_FORTRAN );
    EXPECT( conn( 4, 1 ) == 9 IN_FORTRAN );
    EXPECT( conn( 4, 2 ) == 4 IN_FORTRAN );

    // insert 3 rows with 1 column
    conn.insert( 4, 3, 1 );

    EXPECT( conn.rows() == 9 );
    EXPECT( conn.cols( 4 ) == 1 );
    EXPECT( conn.cols( 5 ) == 1 );
    EXPECT( conn.mincols() == 1 );
    EXPECT( conn.maxcols() == 4 );

    EXPECT( conn( 7, 0 ) == 1 IN_FORTRAN );
    EXPECT( conn( 7, 1 ) == 9 IN_FORTRAN );
    EXPECT( conn( 7, 2 ) == 4 IN_FORTRAN );

    constexpr idx_t cols[3] = {3, 7, 1};
    EXPECT( conn.cols( 2 ) == 4 );
    // insert in position 2, 3 rows with cols[3] number of columns
    conn.insert( 2, 3, cols );

    EXPECT( conn.mincols() == 1 );
    EXPECT( conn.maxcols() == 7 );

    EXPECT( conn.rows() == 12 );
    EXPECT( conn.cols( 2 ) == 3 );
    EXPECT( conn.cols( 3 ) == 7 );
    EXPECT( conn.cols( 4 ) == 1 );
    EXPECT( conn.cols( 5 ) == 4 );

    EXPECT( conn( 5, 0 ) == 4 );
    EXPECT( conn( 5, 1 ) == 13 );
    EXPECT( conn( 5, 2 ) == 55 );
    EXPECT( conn( 5, 3 ) == 78 );
}

CASE( "test_irregular_insert" ) {
    IrregularConnectivity conn( "mesh" );
    EXPECT( conn.rows() == 0 );
    EXPECT( conn.maxcols() == 0 );

    constexpr idx_t vals[4] = {2, 3, 5, 6};
    conn.insert( 0, 1, 4, vals, false );

    EXPECT( conn.rows() == 1 );
    EXPECT( conn.cols( 0 ) == 4 );
    EXPECT( conn.mincols() == 4 );
    EXPECT( conn.maxcols() == 4 );

    EXPECT( conn( 0, 0 ) == 2 );
    EXPECT( conn( 0, 1 ) == 3 );
    EXPECT( conn( 0, 2 ) == 5 );
    EXPECT( conn( 0, 3 ) == 6 );

    EXPECT( conn.row( 0 )( 0 ) == 2 );
    EXPECT( conn.row( 0 )( 1 ) == 3 );
    EXPECT( conn.row( 0 )( 2 ) == 5 );
    EXPECT( conn.row( 0 )( 3 ) == 6 );
}

CASE( "test_block_connectivity" ) {
    idx_t vals[15] = {3, 7, 1, 4, 5, 6, 4, 56, 8, 4, 1, 3, 76, 4, 3};
    BlockConnectivity conn( 3, 5, vals );
    EXPECT( conn.rows() == 3 );
    EXPECT( conn.cols() == 5 );

    EXPECT( conn( 0, 2 ) == 1 );
    EXPECT( conn( 1, 1 ) == 4 );
    EXPECT( conn( 2, 2 ) == 76 );
}

CASE( "test_block_connectivity_add" ) {
    BlockConnectivity conn;
    idx_t vals2[10] = {2, 3, 9, 34, 356, 86, 3, 24, 84, 45};

    conn.add( 2, 5, vals2 );
    EXPECT( conn.rows() == 2 );
    EXPECT( conn.cols() == 5 );

    EXPECT( conn( 0, 2 ) == 9 );
    EXPECT( conn( 0, 4 ) == 356 );
    EXPECT( conn( 1, 1 ) == 3 );
}
CASE( "test_block_connectivity_empty_add" ) {
    BlockConnectivity conn;
    EXPECT( conn.rows() == 0 );
    EXPECT( conn.cols() == 0 );

    idx_t vals2[12] = {2, 3, 9, 34, 356, 86, 3, 24, 84, 45, 2, 2};

    conn.add( 2, 5, vals2 );
    EXPECT( conn.rows() == 2 );
    EXPECT( conn.cols() == 5 );

    EXPECT( conn( 0, 2 ) == 9 );
    EXPECT( conn( 0, 4 ) == 356 );
    EXPECT( conn( 1, 1 ) == 3 );
}

#if 0
CASE("test_multi_block_connectivity_default") {
  idx_t values[22] = {
      1, 3, 4,
      2, 3, 4,
      4, 5, 6, 7,
      23, 54, 6, 9,
      11, 12, 13, 14,
      17, 18,
      21, 24};
  idx_t displ[7] = {0, 3, 6, 10, 14, 18, 20};
  idx_t counts[7] = {3, 3, 4, 4, 4, 2, 2};

  idx_t block_displ[3] = {0, 2, 5};
  idx_t block_cols[3] = {3, 4, 2};
  MultiBlockConnectivity mbc(values, 7, displ, counts, 3, block_displ, block_cols);

  EXPECT(mbc(0, 2) ==  4  );
  EXPECT(mbc(1, 1) ==  3  );
  EXPECT(mbc(2, 2) ==  6  );
  EXPECT(mbc(3, 3) ==  9  );
  EXPECT(mbc(4, 0) ==  11 );
  EXPECT(mbc(5, 0) ==  17 );
  EXPECT(mbc(6, 1) ==  24 );

  EXPECT(mbc(0, 1, 2) ==  4  );
  EXPECT(mbc(1, 1, 0) ==  23 );
  EXPECT(mbc(2, 1, 1) ==  24 );
}
#endif

CASE( "test_multi_block_connectivity_add" ) {
    MultiBlockConnectivity mbc( "mbc" );
    idx_t vals[6] = {1, 3, 4, 2, 3, 4};

    mbc.add( 2, 3, vals, false );
    EXPECT( mbc( 0, 2 ) == 4 );
    EXPECT( mbc( 1, 1 ) == 3 );

    EXPECT( mbc( 0, 1, 2 ) == 4 );

    idx_t vals2[12]{4, 5, 6, 7, 23, 54, 6, 9, 11, 12, 13, 14};
    mbc.add( 3, 4, vals2, false );

    EXPECT( mbc( 2, 2 ) == 6 );
    EXPECT( mbc( 3, 3 ) == 9 );
    EXPECT( mbc( 4, 0 ) == 11 );

    idx_t vals3[4]{17, 18, 21, 24};
    mbc.add( 2, 2, vals3, false );

    EXPECT( mbc( 5, 0 ) == 17 );
    EXPECT( mbc( 6, 1 ) == 24 );

    EXPECT( mbc( 0, 1, 2 ) == 4 );
    EXPECT( mbc( 1, 1, 0 ) == 23 );
    EXPECT( mbc( 2, 1, 1 ) == 24 );
}

CASE( "test_multi_block_connectivity_add_block" ) {
    Log::info() << "\n\n\ntest_multi_block_connectivity_add_block\n" << std::endl;

    MultiBlockConnectivity mbc( "mbc" );

    EXPECT( mbc.blocks() == 0 );

    {
        BlockConnectivity conn1( 3, 5, {3, 7, 1, 4, 5, 6, 4, 56, 8, 4, 1, 3, 76, 4, 3} );
        EXPECT( conn1.owns() );

        EXPECT( conn1( 0, 2 ) == 1 );
        EXPECT( conn1( 1, 1 ) == 4 );
        EXPECT( conn1( 2, 2 ) == 76 );

        mbc.add( conn1 );
    }
    EXPECT( mbc.blocks() == 1 );

    EXPECT( mbc( 0, 2 ) == 1 );
    EXPECT( mbc( 1, 1 ) == 4 );
    EXPECT( mbc( 2, 2 ) == 76 );

    EXPECT( mbc( 0, 0, 2 ) == 1 );
    EXPECT( mbc( 0, 1, 1 ) == 4 );
    EXPECT( mbc( 0, 2, 2 ) == 76 );

    {
        BlockConnectivity conn2( 3, 2, {31, 71, 61, 41, 11, 31} );
        EXPECT( conn2.owns() );
        EXPECT( conn2( 0, 0 ) == 31 );
        EXPECT( conn2( 1, 1 ) == 41 );
        EXPECT( conn2( 2, 0 ) == 11 );

        mbc.add( conn2 );
    }
    EXPECT( mbc.blocks() == 2 );

    EXPECT( mbc( 0, 2 ) == 1 );
    EXPECT( mbc( 1, 1 ) == 4 );
    EXPECT( mbc( 2, 2 ) == 76 );
    EXPECT( mbc( 3, 0 ) == 31 );
    EXPECT( mbc( 4, 1 ) == 41 );
    EXPECT( mbc( 5, 0 ) == 11 );

    EXPECT( mbc( 0, 0, 2 ) == 1 );
    EXPECT( mbc( 0, 1, 1 ) == 4 );
    EXPECT( mbc( 0, 2, 2 ) == 76 );
    EXPECT( mbc( 1, 0, 0 ) == 31 );
    EXPECT( mbc( 1, 1, 1 ) == 41 );
    EXPECT( mbc( 1, 2, 0 ) == 11 );

    const BlockConnectivity& b0 = mbc.block( 0 );
    EXPECT( b0( 0, 2 ) == 1 );
    EXPECT( b0( 1, 1 ) == 4 );
    EXPECT( b0( 2, 2 ) == 76 );

    const BlockConnectivity& b1 = mbc.block( 1 );
    EXPECT( b1( 0, 0 ) == 31 );
    EXPECT( b1( 1, 1 ) == 41 );
    EXPECT( b1( 2, 0 ) == 11 );

    MultiBlockConnectivity& ic = mbc;
    EXPECT( ic( 0, 2 ) == 1 );
    EXPECT( ic( 1, 1 ) == 4 );
    EXPECT( ic( 2, 2 ) == 76 );
    EXPECT( ic( 3, 0 ) == 31 );
    EXPECT( ic( 4, 1 ) == 41 );
    EXPECT( ic( 5, 0 ) == 11 );

    ic.set( 5, 0, 12 );
    EXPECT( b1( 2, 0 ) == 12 );
}

CASE( "test_multi_block_connectivity_add_block_em" ) {
    MultiBlockConnectivity mbc( "mbc" );
    EXPECT( mbc.blocks() == 0 );

    {
        idx_t vals[15]{3, 7, 1, 4, 5, 6, 4, 56, 8, 4, 1, 3, 76, 4, 3};

        BlockConnectivity conn( 3, 5, vals );

        mbc.add( conn );
    }
    EXPECT( mbc.blocks() == 1 );

    EXPECT( mbc( 0, 2 ) == 1 );
    EXPECT( mbc( 1, 1 ) == 4 );
    EXPECT( mbc( 2, 2 ) == 76 );

    EXPECT( mbc( 0, 0, 2 ) == 1 );
    EXPECT( mbc( 0, 1, 1 ) == 4 );
    EXPECT( mbc( 0, 2, 2 ) == 76 );

    {
        idx_t vals5[6]{4, 75, 65, 45, 51, 35};

        BlockConnectivity conn2( 3, 2, vals5 );

        mbc.add( conn2 );
    }
    EXPECT( mbc.blocks() == 2 );

    EXPECT( mbc( 3, 1 ) == 75 );
    EXPECT( mbc( 4, 1 ) == 45 );
    EXPECT( mbc( 5, 0 ) == 51 );

    EXPECT( mbc( 1, 0, 1 ) == 75 );
    EXPECT( mbc( 1, 1, 1 ) == 45 );
    EXPECT( mbc( 1, 2, 0 ) == 51 );
}

CASE( "test_multi_block_connectivity_insert" ) {
    MultiBlockConnectivity mbc( "mbc" );
    EXPECT( mbc.blocks() == 0 );

    mbc.add( BlockConnectivity( 3, 5, {3, 7, 1, 4, 5, 6, 4, 56, 8, 4, 1, 3, 76, 4, 3} ) );

    EXPECT( mbc( 0, 2 ) == 1 );
    EXPECT( mbc( 1, 1 ) == 4 );
    EXPECT( mbc( 2, 2 ) == 76 );

    EXPECT( mbc( 0, 0, 2 ) == 1 );
    EXPECT( mbc( 0, 1, 1 ) == 4 );
    EXPECT( mbc( 0, 2, 2 ) == 76 );

    EXPECT( mbc.block( 0 ).rows() == 3 );

    {
        idx_t vals[10]{31, 71, 61, 41, 42, 11, 31, 33, 54, 56};

        mbc.insert( 0, 2, 5, vals );
    }

    /*
  idx_t check[] = {
      31,71,61,41,42,
      11,31,33,54,56,
      3,7,1,4,5,
      6,4,56,8,4,
      1,3,76,4,3
    };
*/

    EXPECT( mbc.block( 0 ).rows() == 5 );

    EXPECT( mbc( 0, 2 ) == 61 );
    EXPECT( mbc( 1, 1 ) == 31 );

    EXPECT( mbc( 2, 2 ) == 1 );
    EXPECT( mbc( 3, 1 ) == 4 );
    EXPECT( mbc( 4, 2 ) == 76 );

    EXPECT( mbc( 0, 2, 2 ) == 1 );
    EXPECT( mbc( 0, 3, 1 ) == 4 );
    EXPECT( mbc( 0, 4, 2 ) == 76 );

    mbc.add( BlockConnectivity( 3, 2, {4, 75, 65, 45, 51, 35} ) );

    /*
  idx_t check[] = {
    31,71,61,41,42,
    11,31,33,54,56,
    3,7,1,4,5,
    6,4,56,8,4,
    1,3,76,4,3,
    4,75,
    65,45,
    51,35};
*/

    EXPECT( mbc.block( 0 ).rows() == 5 );
    EXPECT( mbc.block( 1 ).rows() == 3 );

    EXPECT( mbc( 0, 2, 2 ) == 1 );
    EXPECT( mbc( 0, 3, 1 ) == 4 );
    EXPECT( mbc( 0, 4, 2 ) == 76 );

    EXPECT( mbc( 5, 1 ) == 75 );
    EXPECT( mbc( 7, 0 ) == 51 );

    EXPECT( mbc( 1, 0, 1 ) == 75 );
    EXPECT( mbc( 1, 2, 0 ) == 51 );

    {
        idx_t vals[10]{19, 72, 62, 42, 43, 12, 32, 34, 55, 57};

        mbc.insert( 5, 2, 5, vals );
    }

    EXPECT( mbc.rows() == 10 );
    EXPECT( mbc.block( 0 ).rows() == 7 );
    EXPECT( mbc.block( 1 ).rows() == 3 );
    EXPECT( mbc.blocks() == 2 );

    /*
  idx_t check[] = {
    31,71,61,41,42,
    11,31,33,54,56,
    3,7,1,4,5,
    6,4,56,8,4,
    1,3,76,4,3,
    19,72,62,42,43, // inserted
    12,32,34,55,57, // inserted

    4,75,
    65,45,
    51,35};
*/

    EXPECT( mbc.block( 0 )( 0, 0 ) == 31 );
    EXPECT( mbc( 0, 0, 0 ) == 31 );
    EXPECT( mbc( 0, 0 ) == 31 );
    EXPECT( mbc.block( 0 )( 5, 0 ) == 19 );
    EXPECT( mbc( 0, 5, 0 ) == 19 );
    EXPECT( mbc( 5, 0 ) == 19 );
    EXPECT( mbc.block( 1 )( 0, 0 ) == 4 );
    EXPECT( mbc( 1, 0, 0 ) == 4 );
    EXPECT( mbc( 7, 0 ) == 4 );

    {
        idx_t vals[] = {67, 38, 89, 2, 5, 7};

        EXPECT_THROWS_AS( mbc.insert( 2, 2, 3, vals ), eckit::AssertionFailed );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
