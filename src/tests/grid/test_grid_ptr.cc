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

#include "atlas/grid/Iterator.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

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

CASE( "test_iterator" ) {
    Grid grid( "O4" );

    for ( atlas::PointXY xy : grid.xy() ) {
        Log::debug() << xy << '\n';
    }

    for ( atlas::PointLonLat ll : grid.lonlat() ) {
        Log::debug() << ll << '\n';
    }
    Log::debug() << std::flush;
}

CASE( "test_iterator_with_predicate" ) {
    Grid grid( "O4" );

    auto filter = []( long n ) {
        bool b = ( n >= 5 && n < 10 );
        if ( b ) {
            Log::debug() << "n = " << n << std::endl;
        }
        return b;
    };

    // auto pred = std::bind( filter, std::placeholders::_1, w_begin, w_end );

    int i( 0 );
    for ( atlas::PointXY xy : grid.xy( filter ) ) {
        Log::debug() << i++ << "  " << xy << std::endl;
    }
    EXPECT( i == 5 );

    Log::debug() << std::flush;
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
