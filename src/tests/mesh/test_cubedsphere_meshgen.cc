/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Tiles.h"
#include "atlas/grid/detail/partitioner/CubedSpherePartitioner.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

CASE( "cubedsphere_mesh_test" ) {
    // Set grid.
    const auto grid = atlas::Grid( "CS-LFR-C-16" );

    // Set mesh config.
    auto meshConfig = util::Config( "partitioner", "equal_regions" );

    // Set mesh generators.
    const auto csMeshgen = atlas::MeshGenerator( "cubedsphere" );              // defaults to cubed sphere partitioner.
    const auto erMeshgen = atlas::MeshGenerator( "cubedsphere", meshConfig );  // Equal regions partitioner.

    const auto csMesh = csMeshgen.generate( grid );
    const auto erMesh = erMeshgen.generate( grid );


    // Set gmsh config.
    auto gmshConfigXy     = atlas::util::Config( "coordinates", "xy" );
    auto gmshConfigXyz    = atlas::util::Config( "coordinates", "xyz" );
    auto gmshConfigLonLat = atlas::util::Config( "coordinates", "lonlat" );

    gmshConfigXy.set( "ghost", true );
    gmshConfigXy.set( "info", true );

    gmshConfigXyz.set( "ghost", true );
    gmshConfigXyz.set( "info", true );

    gmshConfigLonLat.set( "ghost", true );
    gmshConfigLonLat.set( "info", true );

    // Set gmsh objects.
    auto gmshXy     = atlas::output::Gmsh( "cs_xy_mesh.msh", gmshConfigXy );
    auto gmshXyz    = atlas::output::Gmsh( "cs_xyz_mesh.msh", gmshConfigXyz );
    auto gmshLonLat = atlas::output::Gmsh( "cs_lonlat_mesh.msh", gmshConfigLonLat );

    // Write gmsh.
    gmshXy.write( csMesh );
    gmshXyz.write( csMesh );
    gmshLonLat.write( csMesh );

    // Set gmsh objects.
    gmshXy     = atlas::output::Gmsh( "er_xy_mesh.msh", gmshConfigXy );
    gmshXyz    = atlas::output::Gmsh( "er_xyz_mesh.msh", gmshConfigXyz );
    gmshLonLat = atlas::output::Gmsh( "er_lonlat_mesh.msh", gmshConfigLonLat );

    // Write gmsh.
    gmshXy.write( erMesh );
    gmshXyz.write( erMesh );
    gmshLonLat.write( erMesh );
}


}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
