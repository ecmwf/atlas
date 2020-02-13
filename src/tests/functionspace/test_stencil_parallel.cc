/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/parallel/mpi/mpi.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::functionspace;
using namespace atlas::util;
using namespace atlas::grid;


namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test horizontal stencil" ) {
    std::string gridname = eckit::Resource<std::string>( "--grid", "O16" );

    StructuredGrid grid( gridname );
    int halo = eckit::Resource<int>( "--halo", 2 );
    Config config;
    config.set( "halo", halo );
    config.set( "periodic_points", true );
    StructuredColumns fs( grid, Partitioner( "equal_regions" ), config );

    auto gidx = array::make_view<gidx_t, 1>( fs.global_index() );

    HorizontalStencil<4> stencil;

    ComputeHorizontalStencil compute_stencil( grid, stencil.width() );

    auto test = [&]( std::vector<PointXY>& points ) {
        for ( auto p : points ) {
            std::cout << p << std::endl;
            compute_stencil( p.x(), p.y(), stencil );
            for ( idx_t j = 0; j < stencil.width(); ++j ) {
                std::cout << gidx( fs.index( stencil.i( 0, j ), stencil.j( j ) ) ) << std::endl;
                for ( idx_t i = 0; i < stencil.width(); ++i ) {
                    std::cout << stencil.i( i, j ) << ", " << stencil.j( j );
                    std::cout << "   --   "
                              << "x,y = " << fs.compute_xy( stencil.i( i, j ), stencil.j( j ) ) << std::endl;
                    EXPECT( stencil.j( j ) >= fs.j_begin_halo() );
                    EXPECT( stencil.j( j ) < fs.j_end_halo() );
                    EXPECT( stencil.i( i, j ) >= fs.i_begin_halo( stencil.j( j ) ) );
                    EXPECT( stencil.i( i, j ) < fs.i_end_halo( stencil.j( j ) ) );
                }
            }
            std::cout << std::endl;
        }
    };

    std::vector<PointXY> departure_points;

    if ( mpi::comm().size() == 4 ) {
        SECTION( "mpi-size = 4, mpi-rank = 2 " ) {
            if ( mpi::comm().rank() == 2 ) {
                departure_points = {
                    PointXY( 236.25, -30.9375 ),
                };
                test( departure_points );
            }
        }
    }

    if ( mpi::comm().size() == 16 ) {
        SECTION( "ATLAS-186 workaround" ) {
            if ( mpi::comm().rank() == 12 ) {
                // ATLAS-186: point that is found using matching-mesh domain decomposition, which does not really match the halo of StructuredColumns
                departure_points = {
                    PointXY( 205.3125, -67.5 ),
                };
                test( departure_points );
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
