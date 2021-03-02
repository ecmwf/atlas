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

//-----------------------------------------------------------------------------

CASE( "biperiodic_latlon" ) {
    auto& comm       = atlas::mpi::comm();
    const int myproc = comm.rank();


    using namespace atlas::util;
    using namespace atlas;
    using namespace atlas::grid;

    const int Nx = 80, Ny = 80;
    const double xmin = +20, xmax = +60, ymin = +20, ymax = +60;

    std::vector<Spacing> spacings( Ny );

    for ( int i = 0; i < Ny; i++ ) {
        spacings[i] =
            Spacing( Config( "type", "linear" ) | Config( "N", Nx ) | Config( "start", xmin ) | Config( "end", xmax ) );
    }
    StructuredGrid::XSpace xspace( spacings );
    StructuredGrid::YSpace yspace( Config( "type", "linear" ) | Config( "N", Ny ) | Config( "start", ymin ) |
                                   Config( "end", ymax ) );
    Projection proj( Config( "type", "lonlat" ) );

    atlas::StructuredGrid grid( xspace, yspace, proj, Domain() );

    atlas::grid::Distribution dist( grid, atlas::util::Config( "type", "checkerboard" ) );
    atlas::functionspace::StructuredColumns fs( grid, dist,
                                                atlas::util::Config( "halo", 3 ) |
                                                    atlas::util::Config( "periodic_x", true ) |
                                                    atlas::util::Config( "periodic_y", true ) );


    auto f = atlas::Field( "f", atlas::array::DataType::kind<double>(), atlas::array::make_shape( fs.size() ) );

    auto v = atlas::array::make_view<double, 1>( f );

    for ( int i = 0; i < f.size(); i++ ) {
        v( i ) = static_cast<double>( myproc );
    }

    fs.haloExchange( f );

    auto clamp = []( int i, int n ) {
        while ( i < 0 ) {
            i += n;
        }
        while ( i >= n ) {
            i -= n;
        }
        return i;
    };

    auto gv = atlas::array::make_view<atlas::gidx_t, 1>( fs.global_index() );
    auto pv = atlas::array::make_view<int, 1>( fs.partition() );

    for ( int j = fs.j_begin_halo(); j < fs.j_end_halo(); j++ ) {
        for ( int i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); i++ ) {
            int k  = fs.index( i, j );
            int jj = clamp( j, grid.ny() );
            int ii = clamp( i, grid.nx( jj ) );

            int g = grid.index( ii, jj );
            int p = dist.partition( g );

            EXPECT_EQ( v( k ), p );
            EXPECT_EQ( p, pv( k ) );
            EXPECT_EQ( gv( k ) - 1, g );
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
