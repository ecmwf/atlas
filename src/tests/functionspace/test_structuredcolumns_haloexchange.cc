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
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"


#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>


#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

atlas::FieldSet getIJ( const atlas::functionspace::StructuredColumns& fs ) {
    atlas::FieldSet ij;

    auto vi0 = atlas::array::make_view<int, 1>( fs.index_i() );
    auto vj0 = atlas::array::make_view<int, 1>( fs.index_j() );

    auto fi = fs.createField<int>();
    auto fj = fs.createField<int>();

    auto vi1 = atlas::array::make_view<int, 1>( fi );
    auto vj1 = atlas::array::make_view<int, 1>( fj );

    for ( int i = 0; i < fs.size(); i++ ) {
        vi1( i ) = vi0( i );
        vj1( i ) = vj0( i );
    }

    ij.add( fi );
    ij.add( fj );

    return ij;
}


//-----------------------------------------------------------------------------

CASE( "Two haloexchanges for StructuredColumns differing only by 'periodic_points'" ) {
    Grid grid( "L400x200" );

    grid::Distribution dist( grid, grid::Partitioner( "checkerboard" ) );

    functionspace::StructuredColumns fs1( grid, dist, Config( "halo", 1 ) | Config( "periodic_points", true ) );
    functionspace::StructuredColumns fs2( grid, dist, Config( "halo", 1 ) );

    if ( mpi::size() == 1 ) {
        EXPECT_EQ( fs1.sizeOwned(), 80000 );
        EXPECT_EQ( fs2.sizeOwned(), 80000 );
        EXPECT_EQ( fs1.sizeHalo(), 81406 );
        EXPECT_EQ( fs2.sizeHalo(), 81204 );
    }
    if ( mpi::size() == 4 ) {
        EXPECT_EQ( fs1.sizeOwned(), 20000 );
        EXPECT_EQ( fs2.sizeOwned(), 20000 );
        EXPECT_EQ( fs1.sizeHalo(), ( std::vector<int>{20604, 20604, 20604, 20806}[mpi::rank()] ) );
        EXPECT_EQ( fs2.sizeHalo(), ( std::vector<int>{20604, 20604, 20604, 20604}[mpi::rank()] ) );
    }

    auto ij1 = getIJ( fs1 );
    auto ij2 = getIJ( fs2 );

    EXPECT_NO_THROW( fs1.haloExchange( ij1 ) );
    EXPECT_NO_THROW( fs2.haloExchange( ij2 ) );
}

//-----------------------------------------------------------------------------

CASE( "Two haloexchanges for StructuredColumns differing only by 'distribution'" ) {
    Grid grid( "L400x201" );

    grid::Distribution dist1( grid, grid::Partitioner( "checkerboard" ) );
    grid::Distribution dist2( grid, grid::Partitioner( "regular_bands" ) );

    functionspace::StructuredColumns fs1( grid, dist1, Config( "halo", 1 ) );
    functionspace::StructuredColumns fs2( grid, dist2, Config( "halo", 1 ) );

    if ( mpi::size() == 1 ) {
        EXPECT_EQ( fs1.sizeOwned(), 80400 );
        EXPECT_EQ( fs2.sizeOwned(), 80400 );
        EXPECT_EQ( fs1.sizeHalo(), 81606 );
        EXPECT_EQ( fs2.sizeHalo(), 81606 );
    }
    if ( mpi::size() == 4 ) {
        EXPECT_EQ( fs1.sizeOwned(), 20100 );
        EXPECT_EQ( fs2.sizeOwned(), ( std::vector<int>{20400, 20000, 20000, 20000}[mpi::rank()] ) );
        EXPECT_EQ( fs1.sizeHalo(), ( std::vector<int>{20706, 20706, 20706, 20706}[mpi::rank()] ) );
        EXPECT_EQ( fs2.sizeHalo(), ( std::vector<int>{21306, 20904, 20904, 20904}[mpi::rank()] ) );
    }

    auto ij1 = getIJ( fs1 );
    auto ij2 = getIJ( fs2 );

    EXPECT_NO_THROW( fs1.haloExchange( ij1 ) );
    EXPECT_NO_THROW( fs2.haloExchange( ij2 ) );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
