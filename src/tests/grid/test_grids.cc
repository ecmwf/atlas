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

#include "atlas/grid.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------


CASE( "test_ij2gidx" ) {
    StructuredGrid n16 = Grid( "N16" );

    for ( int j = 0, jglo = 0; j < n16.ny(); j++ ) {
        for ( int i = 0; i < n16.nx( j ); i++, jglo++ ) {
            idx_t i1, j1;
            n16.index2ij( jglo, i1, j1 );
            EXPECT( n16.index( i, j ) == jglo );
            EXPECT( i1 == i );
            EXPECT( j1 == j );
        }
    }
}

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

CASE( "ATLAS-255: regular Gaussian grid with global domain" ) {
    GlobalDomain globe;
    Grid grid( "F80", globe );
    EXPECT( GaussianGrid( grid ) );
}

CASE( "ATLAS-255: reduced Gaussian grid with global domain" ) {
    GlobalDomain globe;
    Grid grid( "O80", globe );
    EXPECT( GaussianGrid( grid ) );
}

CASE( "test_from_string_L32" ) {
    Grid grid;
    EXPECT( not grid );

    grid = Grid( "L32" );
    EXPECT( grid );
    EXPECT( StructuredGrid( grid ) == true );
    EXPECT( RegularGrid( grid ) == true );

    auto structured = StructuredGrid( grid );
    EXPECT( structured.ny() == 65 );
    EXPECT( structured.periodic() == true );
    EXPECT( structured.nx( 0 ) == 128 );
    EXPECT( structured.y().front() == 90. );
    EXPECT( structured.y().back() == -90. );

    auto regular = RegularGrid( grid );
    EXPECT( regular.ny() == 65 );
    EXPECT( regular.periodic() == true );
    EXPECT( regular.nx() == 128 );
    EXPECT( regular.y().front() == 90. );
    EXPECT( regular.y().back() == -90. );
}

CASE( "test_from_string_O32" ) {
    Grid grid;
    EXPECT( not grid );

    grid = Grid( "O32" );
    EXPECT( grid );

    EXPECT( StructuredGrid( grid ) );
    EXPECT( !RegularGrid( grid ) );

    auto structured = StructuredGrid( grid );
    EXPECT( structured.ny() == 64 );
    EXPECT( structured.periodic() == true );
    EXPECT( structured.nx().front() == 20 );
}

CASE( "test_from_string_O32_with_domain" ) {
    Grid grid;
    EXPECT( not grid );

    grid = Grid( "O32", RectangularDomain( {0, 90}, {0, 90} ) );
    EXPECT( grid );

    EXPECT( StructuredGrid( grid ) );
    EXPECT( !RegularGrid( grid ) );

    auto structured = StructuredGrid( grid );
    EXPECT( structured.ny() == 32 );
    EXPECT( structured.periodic() == false );
    EXPECT( structured.nx().front() == 6 );

    output::Gmsh gmsh( "test_grid_ptr_O32_subdomain.msh" );
    Mesh mesh = StructuredMeshGenerator().generate( grid );
    gmsh.write( mesh );
}

CASE( "test_structured_1" ) {
    std::stringstream json;
    json << "{"
            "\"type\" : \"structured\","
            "\"yspace\" : { \"type\":\"linear\", \"N\":9,  \"start\":90, "
            "\"end\":-90 },"
            "\"xspace\" : { \"type\":\"linear\", \"N\":16, \"start\":0,  "
            "\"end\":360, \"endpoint\":false }"
            "}";
    json.seekp( 0 );

    Grid grid;
    EXPECT( not grid );

    Config json_config( json );

    grid = StructuredGrid( json_config );
    EXPECT( grid );
    EXPECT( StructuredGrid( grid ) );
    EXPECT( RegularGrid( grid ) );

    auto structured = StructuredGrid( grid );
    EXPECT( structured.ny() == 9 );
    EXPECT( structured.periodic() == true );
    EXPECT( structured.nx( 0 ) == 16 );
    EXPECT( structured.y().front() == 90. );
    EXPECT( structured.y().back() == -90. );

    auto regular = RegularGrid( grid );
    EXPECT( regular.ny() == 9 );
    EXPECT( regular.periodic() == true );
    EXPECT( regular.nx() == 16 );
    EXPECT( regular.y().front() == 90. );
    EXPECT( regular.y().back() == -90. );

    output::Gmsh gmsh( "test_grid_ptr.msh" );
    Mesh mesh = StructuredMeshGenerator().generate( grid );
    gmsh.write( mesh );
}

CASE( "test_structured_2" ) {
    using XSpace     = StructuredGrid::XSpace;
    using YSpace     = StructuredGrid::YSpace;
    using Domain     = StructuredGrid::Domain;
    using Projection = StructuredGrid::Projection;
    StructuredGrid grid( XSpace( {0., 360.}, {2, 4, 6, 6, 4, 2}, false ),
                         YSpace( grid::LinearSpacing( {90., -90.}, 6 ) ), Projection(), Domain() );
    EXPECT( grid );

    output::Gmsh gmsh( "test_grid_ptr_structured_2.msh" );
    Mesh mesh = StructuredMeshGenerator().generate( grid );
    gmsh.write( mesh );

    Log::info() << grid.spec() << std::endl;

    Grid newgrid( grid.spec() );
    Log::info() << newgrid.spec() << std::endl;

    Log::info() << "original: " << grid.uid() << std::endl;
    Log::info() << "fromspec: " << newgrid.uid() << std::endl;
    EXPECT( grid == newgrid );
}

CASE( "test_structured_3" ) {
    StructuredGrid grid( "O32" );
    EXPECT( grid );

    Log::info() << grid.spec() << std::endl;

    Grid newgrid( grid.spec() );
    Log::info() << newgrid.spec() << std::endl;

    Log::info() << "original: " << grid.uid() << std::endl;
    Log::info() << "fromspec: " << newgrid.uid() << std::endl;
    EXPECT( grid == newgrid );
    EXPECT( grid.name() == "O32" );
    EXPECT( newgrid.name() == "O32" );
}

CASE( "test_structured_triangulated" ) {
    Grid grid;

    // Create grid
    {
        using XSpace = StructuredGrid::XSpace;
        using YSpace = StructuredGrid::YSpace;
        auto xspace  = util::Config{};
        xspace.set( "type", "linear" );
        xspace.set( "N", 16 );
        xspace.set( "length", 360 );
        xspace.set( "endpoint", false );
        xspace.set( "start[]", []() {
            auto startpts = std::vector<double>( 8 );
            for ( int i = 0; i < 8; ++i ) {
                startpts[i] = i * 12.;
            }
            return startpts;
        }() );
        grid = StructuredGrid{XSpace{xspace}, YSpace{grid::LinearSpacing{{90., -90.}, 8}}};
    }

    EXPECT( grid );
    Log::info() << grid.spec() << std::endl;

    // Create and output mesh
    {
        auto meshgen = StructuredMeshGenerator{util::Config( "angle", -1. )};
        auto mesh    = meshgen.generate( grid );
        auto gmsh    = output::Gmsh{"structured_triangulated.msh"};
        gmsh.write( mesh );
    }
}

CASE( "test_structured_from_config" ) {
    Config config;
    config.set( "type", "structured" );
    config.set( "xspace", []() {
        Config config;
        config.set( "type", "linear" );
        config.set( "N", 40 );
        config.set( "start", 5 );
        config.set( "end", 365 );
        config.set( "endpoint", false );
        return config;
    }() );
    config.set( "yspace", []() {
        Config config;
        config.set( "type", "custom" );
        config.set( "N", 9 );
        config.set( "values", std::vector<double>{5., 15., 25., 35., 45., 55., 65., 75., 85.} );
        return config;
    }() );
    StructuredGrid g{config};
    for ( idx_t j = 0; j < g.ny(); ++j ) {
        EXPECT_EQ( g.nx( j ), 40 );
        EXPECT_EQ( g.x( 0, j ), 5. );
        EXPECT_EQ( g.x( g.nx( j ), j ), 365. );
        EXPECT_EQ( g.dx( j ), 9. );
        EXPECT_EQ( g.xmin( j ), 5. );
    }
    EXPECT( not g.domain().global() );
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
