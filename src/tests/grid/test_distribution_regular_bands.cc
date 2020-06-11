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

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

atlas::FieldSet getIJ( const atlas::functionspace::StructuredColumns& fs ) {
    atlas::FieldSet ij;

    //    auto i = atlas::array::make_view<int, 1>( fs.index_i() );
    //    auto j = atlas::array::make_view<int, 1>( fs.index_j() );

    ij.add( fs.index_i() );
    ij.add( fs.index_j() );

    return ij;
}


//-----------------------------------------------------------------------------

CASE( "test_bands" ) {
    int nproc = mpi::size();

    StructuredGrid grid = Grid( "L200x101" );
    grid::Distribution dist( grid, atlas::util::Config( "type", "checkerboard" ) | Config( "bands", nproc ) );

    for ( int i = 1; i < grid.size(); i++ ) {
        EXPECT( dist.partition( i - 1 ) <= dist.partition( i ) );
    }
}

CASE( "test_regular_bands" ) {
    int nproc = mpi::size();

    std::vector<std::string> gridnames = {"L40x21", "L40x20", "Slat100x50"};

    for ( auto gridname : gridnames ) {
        SECTION( gridname ) {
            auto grid    = RegularGrid( gridname );
            const int nx = grid.nx();
            const int ny = grid.ny();

            bool equivalent_with_checkerboard =
                ( ny % nproc ) == 0;  // Compare regular band & checkerboard distributions when possible

            grid::Distribution checkerboard( grid, grid::Partitioner( "checkerboard", Config( "bands", nproc ) ) );
            grid::Distribution regularbands( grid, grid::Partitioner( "regular_bands" ) );

            EXPECT( regularbands.footprint() < 100 );

            const auto& nb_pts2 = regularbands.nb_pts();

            int count = 0;
            for ( int i = 0; i < grid.size(); i++ ) {
                if ( regularbands.partition( i ) == mpi::rank() ) {
                    count++;
                }
            }

            EXPECT_EQ( count, nb_pts2[mpi::rank()] );


            if ( equivalent_with_checkerboard ) {
                for ( int i = 0; i < grid.size(); i++ ) {
                    EXPECT_EQ( checkerboard.partition( i ), regularbands.partition( i ) );
                }
            }
            else {
                Log::warning() << "WARNING: checkerboard not expected to be equal!" << std::endl;
            }

            // Expect each points with constant latitude to be on the same partition as the first on each latitude.
            for ( int iy = 0, jglo = 0; iy < ny; iy++ ) {
                int jglo0 = jglo;
                for ( int ix = 0; ix < nx; ix++, jglo++ ) {
                    EXPECT_EQ( regularbands.partition( jglo ), regularbands.partition( jglo0 ) );
                }
            }

            functionspace::StructuredColumns fs_checkerboard( grid, checkerboard,
                                                              Config( "halo", 1 ) | Config( "periodic_points", true ) );
            functionspace::StructuredColumns fs_regularbands( grid, regularbands,
                                                              Config( "halo", 1 ) | Config( "periodic_points", true ) );

            auto ij1 = getIJ( fs_checkerboard );
            auto ij2 = getIJ( fs_regularbands );

            fs_checkerboard.haloExchange( ij1 );
            fs_regularbands.haloExchange( ij2 );

            if ( equivalent_with_checkerboard ) {
                EXPECT_EQ( fs_regularbands.size(), fs_checkerboard.size() );
                EXPECT_EQ( fs_regularbands.sizeOwned(), fs_checkerboard.sizeOwned() );

                auto i1 = array::make_view<int, 1>( ij1[0] );
                auto j1 = array::make_view<int, 1>( ij1[1] );
                auto i2 = array::make_view<int, 1>( ij2[0] );
                auto j2 = array::make_view<int, 1>( ij2[1] );

                for ( int k = 0; k < fs_checkerboard.sizeOwned(); k++ ) {
                    EXPECT_EQ( i1[k], i2[k] );
                    EXPECT_EQ( j1[k], j2[k] );
                }
            }

            for ( int j = fs_regularbands.j_begin_halo(); j < fs_regularbands.j_end_halo(); j++ ) {
                EXPECT_EQ( fs_regularbands.i_begin_halo( j ), -1 );
                EXPECT_EQ( fs_regularbands.i_end_halo( j ), nx + 2 );
            }
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
