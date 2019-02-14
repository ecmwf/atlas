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

#include "eckit/types/FloatCompare.h"

#include "atlas/grid.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using Grid = atlas::Grid;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_factory" ) {
    StructuredGrid structured = Grid( "N80" );

    Grid grid = Grid( "N24" );

    std::cout << "structured.ny() = " << structured.ny() << std::endl;
    std::cout << "grid.npts() = " << grid.size() << std::endl;
}

CASE( "test_regular_gg" ) {
    RegularGrid grid( "F32" );

    EXPECT( grid.ny() == 64 );
    EXPECT( grid.size() == 8192 );
    // EXPECT(grid.type() == "regular_gaussian");

    // Full Gaussian Grid

    Grid::Config config;
    config.set( "type", "regular_gaussian" );
    config.set( "N", 32 );
    grid = Grid( config );
    EXPECT( grid.size() == 8192 );
    // EXPECT(grid.type() == "regular_gaussian");
}

CASE( "test_reduced_gg" ) {
    StructuredGrid grid;

    grid = Grid( "N32" );
    EXPECT( grid.ny() == 64 );
    EXPECT( grid.size() == 6114 );

    grid = ReducedGaussianGrid( {4, 6, 8, 8, 6, 4} );

    EXPECT( grid.ny() == 6 );
    EXPECT( grid.size() == 8 + 12 + 16 );
}

CASE( "test_reduced_gg_ifs" ) {
    StructuredGrid grid( "N32" );

    // EXPECT(grid.N() ==    32);
    EXPECT( grid.ny() == 64 );
    EXPECT( grid.size() == 6114 );
    // EXPECT(grid.type() == "classic_gaussian");
}

CASE( "test_regular_ll" ) {
    // Constructor for N=8
    idx_t nlon = 32;
    idx_t nlat = 16;
    std::stringstream name;
    name << "Slat" << nlon << "x" << nlat;
    RegularGrid grid( name.str() );

    EXPECT( grid.nx() == nlon );
    EXPECT( grid.ny() == nlat );
    EXPECT( grid.size() == 512 );
    // EXPECT(grid.type() == "shifted_lat");
    EXPECT( is_approximately_equal( grid.y( 0 ), 90. - 0.5 * ( 180. / 16. ) ) );
    EXPECT( is_approximately_equal( grid.y( grid.ny() - 1 ), -90. + 0.5 * ( 180. / 16. ) ) );
    EXPECT( is_approximately_equal( grid.x( 0 ), 0. ) );
    EXPECT( is_approximately_equal( grid.x( grid.nx() - 1 ), 360. - 360. / 32. ) );

    // Construct using builders/factories

    // Global Grid
    Grid::Config config1;
    config1.set( "type", "shifted_lat" );
    config1.set( "nx", 32 );
    config1.set( "ny", 16 );
    grid = Grid( config1 );
    EXPECT( grid.size() == 512 );
    // EXPECT(gridptr->type() == "shifted_lat");

    Grid::Config config2;
    config2.set( "type", "shifted_lat" );
    config2.set( "N", 8 );
    grid = Grid( config2 );
    EXPECT( grid.size() == 512 );
    // EXPECT(gridptr->type() == "shifted_lat");

    RegularGrid ll_poles( "L4x3" );
    EXPECT( ll_poles.nx() == 4 );
    EXPECT( ll_poles.ny() == 3 );

    RegularGrid ll_nopoles( "Slat4x2" );
    EXPECT( ll_nopoles.nx() == 4 );
    EXPECT( ll_nopoles.ny() == 2 );
    EXPECT( is_approximately_equal( ll_nopoles.y( 0 ), 45. ) );   // tolerance was previously 1.e-5
    EXPECT( is_approximately_equal( ll_nopoles.y( 1 ), -45. ) );  // tolerance was previously 1.e-5
    EXPECT( is_approximately_equal( ll_nopoles.x( 0 ), 0. ) );    // tolerance was previously 1.e-5
    EXPECT( is_approximately_equal( ll_nopoles.x( 1 ), 90. ) );   // tolerance was previously 1.e-5
}

CASE( "test_reducedgaussian" ) {
    StructuredGrid N640( "N640" );
    EXPECT( N640.size() == 2140702 );
    ReducedGaussianGrid custom( N640.nx() );
    EXPECT( N640.size() == custom.size() );
}

CASE( "test_cropping previous case" ) {
    StructuredGrid grid( "N32" );
    EXPECT( grid.ny() == 64 );
    EXPECT( grid.size() == 6114 );

    StructuredGrid cropped( grid, RectangularDomain( {-27, 45}, {33, 73} ) );
    EXPECT( cropped.ny() == 14 );
    EXPECT( cropped.size() == 267 );
}

CASE( "cropping with line at north pole" ) {
    StructuredGrid grid( "L16", RectangularDomain( {0, 360}, {90, 90} ) );
    EXPECT( grid.ny() == 1 );
    EXPECT( grid.nx( 0 ) == 64 );
    EXPECT( grid.size() == 64 );
}

CASE( "cropping with line at south pole" ) {
    StructuredGrid grid( "L16", RectangularDomain( {0, 360}, {-90, -90} ) );
    EXPECT( grid.ny() == 1 );
    EXPECT( grid.nx( 0 ) == 64 );
    EXPECT( grid.size() == 64 );
}

CASE( "cropping with line at equator" ) {
    StructuredGrid grid( "L16", RectangularDomain( {0, 360}, {0, 0} ) );
    EXPECT( grid.ny() == 1 );
    EXPECT( grid.nx( 0 ) == 64 );
    EXPECT( grid.size() == 64 );
}

CASE( "cropping single point at equator" ) {
    StructuredGrid grid( "L16", RectangularDomain( {0, 0}, {0, 0} ) );
    EXPECT( grid.ny() == 1 );
    EXPECT( grid.nx( 0 ) == 1 );
    EXPECT( grid.size() == 1 );
}

CASE( "Create cropped unstructured grid using rectangular domain" ) {
    StructuredGrid agrid( "L8" );
    auto domain = RectangularDomain( {-27, 45}, {33, 73} );
    StructuredGrid sgrid( agrid, domain );
    UnstructuredGrid ugrid( agrid, domain );
    EXPECT( ugrid.size() == sgrid.size() );
}

CASE( "Create cropped unstructured grid using zonal domain" ) {
    StructuredGrid agrid( "L8" );
    auto domain = ZonalBandDomain( {33, 73} );
    StructuredGrid sgrid( agrid, domain );
    UnstructuredGrid ugrid( agrid, domain );
    EXPECT( ugrid.size() == sgrid.size() );
}

CASE( "Create unstructured from unstructured" ) {
    StructuredGrid agrid( "L8" );
    UnstructuredGrid global_unstructured( agrid, Domain() );
    EXPECT( UnstructuredGrid( global_unstructured ) );
    EXPECT( global_unstructured.size() == agrid.size() );
    auto domain = ZonalBandDomain( {33, 73} );
    UnstructuredGrid ugrid( global_unstructured, domain );
    EXPECT( ugrid.size() == StructuredGrid( agrid, domain ).size() );
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
