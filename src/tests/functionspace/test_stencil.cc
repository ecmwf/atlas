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
#include "eckit/memory/ScopedPtr.h"
#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/MicroDeg.h"

#include "Stencil.h"
#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------


CASE( "test finding of North-West grid point" ) {
    std::string gridname = eckit::Resource<std::string>( "--grid", "O8" );

    grid::StructuredGrid grid( gridname );

    constexpr double tol = 0.5e-6;

    ComputeNorth compute_j_north( grid );
    ComputeWest compute_i_west( grid );

    struct IJ {
        idx_t i;
        idx_t j;
    };
    idx_t ny = grid.ny();

    if ( mpi::comm().size() == 1 ) {
        auto entries = {
            std::make_tuple( PointXY{0. + 0.5 * tol, grid.y( 0 ) + 0.5 * tol}, IJ{0, 0} ),
            std::make_tuple( PointXY{0. - 0.5 * tol, grid.y( 0 ) - 0.5 * tol}, IJ{0, 0} ),
            std::make_tuple( PointXY{0. + 2.0 * tol, grid.y( 0 ) + 2.0 * tol}, IJ{0, -1} ),
            std::make_tuple( PointXY{0. - 2.0 * tol, grid.y( 0 ) - 2.0 * tol}, IJ{-1, 0} ),
            std::make_tuple( PointXY{360. + 0.5 * tol, grid.y( ny - 1 ) + 0.5 * tol}, IJ{grid.nx( ny - 1 ), ny - 1} ),
            std::make_tuple( PointXY{360. - 0.5 * tol, grid.y( ny - 1 ) - 0.5 * tol}, IJ{grid.nx( ny - 1 ), ny - 1} ),
            std::make_tuple( PointXY{360. + 2.0 * tol, grid.y( ny - 1 ) + 2.0 * tol}, IJ{grid.nx( ny - 2 ), ny - 2} ),
            std::make_tuple( PointXY{360. - 2.0 * tol, grid.y( ny - 1 ) - 2.0 * tol},
                             IJ{grid.nx( ny - 1 ) - 1, ny - 1} ),
        };
        for ( auto entry : entries ) {
            auto p     = std::get<0>( entry );
            auto index = std::get<1>( entry );
            EXPECT( compute_j_north( p.y() ) == index.j );
            EXPECT( compute_i_west( p.x(), index.j ) == index.i );
        }
    }
}

CASE( "test horizontal stencil" ) {
    //if ( mpi::comm().size() == 1 ) {
    std::string gridname = eckit::Resource<std::string>( "--grid", "O8" );

    grid::StructuredGrid grid( gridname );
    int halo = eckit::Resource<int>( "--halo", 2 );
    util::Config config;
    config.set( "halo", halo );
    config.set( "levels", 9 );
    config.set( "periodic_points", true );
    functionspace::StructuredColumns fs( grid, grid::Partitioner( "equal_regions" ), config );

    double tol = 0.5e-6;

    constexpr int stencil_width = 4;
    HorizontalStencil<stencil_width> stencil;

    ComputeHorizontalStencil compute_stencil( grid, stencil.width() );

    auto departure_points = {
        PointXY( 0., 90. ),
        PointXY( 0., -90. ),
        PointXY( 0., 0. ),
        PointXY( 360., 0. ),
    };
    for ( auto p : departure_points ) {
        Log::info() << p << std::endl;
        compute_stencil( p.x(), p.y(), stencil );
        for ( idx_t j = 0; j < stencil.width(); ++j ) {
            Log::info() << stencil.i( j ) << " " << stencil.j( j ) << "   --   "
                        << "x,y = " << fs.compute_xy( stencil.i( j ), stencil.j( j ) ) << std::endl;
        }
        Log::info() << std::endl;
    }
}

CASE( "test vertical stencil" ) {
    idx_t nlev = 10;
    std::vector<double> zcoord( nlev + 2 );
    double dzcoord = 1. / double( nlev + 1 );
    for ( idx_t jlev = 0; jlev <= nlev + 1; ++jlev ) {
        zcoord[jlev] = jlev * dzcoord;
    }
}


//-----------------------------------------------------------------------------


CASE( "ifs method to find nearest grid point" ) {
    // see satrad/module/gaussgrid.F90
    std::string gridname = eckit::Resource<std::string>( "--grid", "O8" );
    grid::StructuredGrid grid( gridname );

    auto p = PointXY{0., grid.y( 0 )};
    idx_t kgrib_lat, kgrib_lon;
    {
        double x = p.x();
        double y = p.y();
        std::vector<double> pdlat( grid.ny() );
        for ( idx_t j = 0; j < grid.ny(); ++j ) {
            pdlat[j] = std::abs( y - grid.y( j ) );
        }

        auto iterator = std::min_element( pdlat.begin(), pdlat.end() );
        kgrib_lat     = ( iterator - pdlat.begin() );

        double zfirstlon = grid.x( 0, kgrib_lat );
        double zdlon     = grid.x( 1, kgrib_lat ) - zfirstlon;
        double zsafelon  = std::fmod( x - zfirstlon + 720., 360. );
        kgrib_lon        = std::fmod( std::round( zsafelon / zdlon ), grid.nx( kgrib_lat ) );
    }
    EXPECT( kgrib_lon == 0 );
    EXPECT( kgrib_lat == 0 );
}

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
