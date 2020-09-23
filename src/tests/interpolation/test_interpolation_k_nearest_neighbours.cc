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
#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::NodeColumns;
using atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_multiple_fs" ) {
    Grid grid1( "L90x45" );
    Grid grid2( "O8" );

    Mesh mesh1 = StructuredMeshGenerator().generate( grid1 );
    Mesh mesh2 = StructuredMeshGenerator().generate( grid2 );

    functionspace::NodeColumns fs11( mesh1, option::halo( 1 ) );
    functionspace::NodeColumns fs12( mesh1, option::halo( 2 ) );

    auto fs1 = fs11;
    functionspace::NodeColumns fs2( mesh2, option::halo( 1 ) );

    Interpolation interpolation12( Config( "type", "k-nearest-neighbours" ) | Config( "k-nearest-neighbours", 5 ), fs1,
                                   fs2 );

    auto f1 = fs1.createField<double>( Config( "name", "source" ) );
    auto f2 = fs2.createField<double>( Config( "name", "target" ) );

    auto v1 = array::make_view<double,1>( f1 );
    v1.assign( 1. );

    interpolation12.execute( f1, f2 );
}

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
