/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


/// Linked JIRA issue: ATLAS-247
/// For now the CellColumns functionspace only works when the mesh has elements of only one type.
/// For this reason we triangulate the mesh always, and don't patch the pole
/// The problem is in the computation of mesh.cells().remote_index() during BuildHalo called
/// within the CellColumns constructor. Meshes without any parallel halo will also succeed as
/// the BuildHalo routine is not called.


#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/mpi/mpi.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

Mesh generate_mesh() {
    auto grid = Grid{"O16"};
    auto meshgenerator =
        MeshGenerator{util::Config( "type", "structured" )( "patch_pole", false )( "triangulate", true )};
    return meshgenerator.generate( grid );
}

void set_field_values( const Mesh& mesh, Field& field ) {
    auto value     = array::make_view<int, 1>( field );
    auto partition = array::make_view<int, 1>( mesh.cells().partition() );
    auto halo      = array::make_view<int, 1>( mesh.cells().halo() );

    EXPECT( field.shape( 0 ) == mesh.cells().size() );

    const size_t nb_cells = mesh.cells().size();
    for ( size_t j = 0; j < nb_cells; ++j ) {
        if ( halo( j ) ) {
            value( j ) = -1;
        }
        else {
            value( j ) = partition( j );
        }
    }
}

void check_field_values( const Mesh& mesh, Field& field ) {
    auto value            = array::make_view<int, 1>( field );
    auto partition        = array::make_view<int, 1>( mesh.cells().partition() );
    const size_t nb_cells = mesh.cells().size();
    for ( size_t j = 0; j < nb_cells; ++j ) {
        EXPECT( value( j ) == partition( j ) );
    }
}

//-----------------------------------------------------------------------------

CASE( "test_functionspace_CellColumns_no_halo" ) {
    Mesh mesh = generate_mesh();
    CellColumns fs( mesh );

    Field field( fs.createField<int>() );

    set_field_values( mesh, field );

    fs.haloExchange( field );

    check_field_values( mesh, field );
}

CASE( "test_functionspace_CellColumns_halo_1" ) {
    Mesh mesh = generate_mesh();
    CellColumns fs( mesh, option::halo( 1 ) );

    Field field( fs.createField<int>() );

    set_field_values( mesh, field );

    fs.haloExchange( field );

    check_field_values( mesh, field );

    output::Gmsh output( "cellcolumns_halo1.msh" );
    output.write( mesh, util::Config( "ghost", true ) );
    output.write( field );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
