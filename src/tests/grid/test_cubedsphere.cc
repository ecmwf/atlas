/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Tiles.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

CASE( "cubedsphere_tile_test" ) {
    auto tileConfig1 = atlas::util::Config( "type", "cubedsphere_lfric" );
    auto lfricTiles  = atlas::grid::CubedSphereTiles( tileConfig1 );
    EXPECT( lfricTiles.type() == "cubedsphere_lfric" );

    auto tileConfig2 = atlas::util::Config( "type", "cubedsphere_fv3" );
    auto fv3Tiles    = atlas::grid::CubedSphereTiles( tileConfig2 );
    EXPECT( fv3Tiles.type() == "cubedsphere_fv3" );

    auto lfricTiles2 = atlas::grid::CubedSphereTiles( "cubedsphere_lfric" );
    EXPECT( lfricTiles2.type() == "cubedsphere_lfric" );

    auto fv3Tiles2 = atlas::grid::CubedSphereTiles( "cubedsphere_fv3" );
    EXPECT( fv3Tiles.type() == "cubedsphere_fv3" );
}


//-----------------------------------------------------------------------------

CASE( "test_iterator" ) {
    std::vector<int> resolutions{1, 2, 4, 8};
    std::vector<std::string> grid_prefixes{"CS-EA-L-", "CS-ED-L-", "CS-EA-C-", "CS-ED-C-", "CS-LFR-L-", "CS-LFR-C-"};


    for ( auto resolution : resolutions ) {
        for ( auto& grid_prefix : grid_prefixes ) {
            std::string grid_name = grid_prefix + std::to_string( resolution );
            SECTION( grid_name ) {
                if ( grid_name == "CS-LFR-L-1" ) {
                    Log::error() << eckit::Colour::red << "TODO: Fix me!!!. Skipping..." << eckit::Colour::reset
                                 << std::endl;
                    continue;
                }
                Grid g( grid_name );

                // Test xy
                {
                    std::vector<PointXY> coordinates_1;
                    std::vector<PointXY> coordinates_2;
                    std::vector<PointXY> coordinates_3;
                    {
                        for ( auto crd : g.xy() ) {
                            coordinates_1.push_back( crd );
                        }
                    }
                    {
                        auto iterator = g.xy().begin();
                        for ( int n = 0; n < g.size(); ++n ) {
                            coordinates_2.push_back( *iterator );
                            iterator += 1;
                        }
                    }
                    {
                        auto iterator = g.xy().begin();
                        PointXY crd;
                        while ( iterator.next( crd ) ) {
                            coordinates_3.push_back( crd );
                        }
                    }
                    EXPECT_EQ( coordinates_1.size(), g.size() );
                    EXPECT_EQ( coordinates_2.size(), g.size() );
                    EXPECT_EQ( coordinates_3.size(), g.size() );
                    EXPECT_EQ( coordinates_2, coordinates_1 );
                    EXPECT_EQ( coordinates_3, coordinates_1 );
                }

                // Test lonlat
                {
                    std::vector<PointLonLat> coordinates_1;
                    std::vector<PointLonLat> coordinates_2;
                    std::vector<PointLonLat> coordinates_3;
                    {
                        for ( auto crd : g.lonlat() ) {
                            coordinates_1.push_back( crd );
                        }
                    }
                    {
                        auto iterator = g.lonlat().begin();
                        for ( int n = 0; n < g.size(); ++n ) {
                            coordinates_2.push_back( *iterator );
                            iterator += 1;
                        }
                    }
                    {
                        auto iterator = g.lonlat().begin();
                        PointLonLat crd;
                        while ( iterator.next( crd ) ) {
                            coordinates_3.push_back( crd );
                        }
                    }
                    EXPECT_EQ( coordinates_1.size(), g.size() );
                    EXPECT_EQ( coordinates_2.size(), g.size() );
                    EXPECT_EQ( coordinates_3.size(), g.size() );
                    EXPECT_EQ( coordinates_2, coordinates_1 );
                    EXPECT_EQ( coordinates_3, coordinates_1 );
                }
            }
        }
    }
}


CASE( "cubedsphere_grid_mesh_field_test" ) {
    // THIS IS TEMPORARY!
    // I expect this will be replaced by some more aggressive tests.

    // Set grid.
    const auto grid = atlas::Grid( "CS-EA-L-2" );

    atlas::Log::info() << grid->type() << std::endl;
    atlas::Log::info() << grid.size() << std::endl;


    // Set mesh.
    auto meshGen = atlas::MeshGenerator( "cubedsphere" );
    auto mesh    = meshGen.generate( grid );

    // Set functionspace
    auto functionSpace = atlas::functionspace::NodeColumns( mesh );

    auto ghostIdx = mesh.nodes().metadata().get<std::vector<idx_t>>( "ghost-global-idx" );
    auto ownedIdx = mesh.nodes().metadata().get<std::vector<idx_t>>( "owned-global-idx" );

    // Print out ghost global indices with corresponding owned global indices
    auto ownedIdxIt = ownedIdx.begin();
    for ( auto iGhost : ghostIdx ) {
        std::cout << iGhost << " " << *ownedIdxIt++ << std::endl;
    }

    // Set field
    auto field = functionSpace.ghost();

    // Set gmsh config.
    auto gmshConfigXy     = atlas::util::Config( "coordinates", "xy" ) | atlas::util::Config( "ghost", false );
    auto gmshConfigXyz    = atlas::util::Config( "coordinates", "xyz" ) | atlas::util::Config( "ghost", false );
    auto gmshConfigLonLat = atlas::util::Config( "coordinates", "lonlat" ) | atlas::util::Config( "ghost", false );

    // Set source gmsh object.
    const auto gmshXy     = atlas::output::Gmsh( "cs_xy_mesh.msh", gmshConfigXy );
    const auto gmshXyz    = atlas::output::Gmsh( "cs_xyz_mesh.msh", gmshConfigXyz );
    const auto gmshLonLat = atlas::output::Gmsh( "cs_lonlat_mesh.msh", gmshConfigLonLat );

    // Write gmsh.
    gmshXy.write( mesh );
    gmshXy.write( field );
    gmshXyz.write( mesh );
    gmshXyz.write( field );
    gmshLonLat.write( mesh );
    gmshLonLat.write( field );
}


}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
