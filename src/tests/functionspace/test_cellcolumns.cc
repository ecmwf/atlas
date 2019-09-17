/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/Library.h"
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
    auto grid          = Grid{"O16"};
    auto meshgenerator = StructuredMeshGenerator{util::Config( "patch_pole", false )( "triangulate", true )};
    return meshgenerator.generate( grid );
}

//-----------------------------------------------------------------------------

CASE( "test_functionspace_CellColumns_no_halo" ) {
    Mesh mesh = generate_mesh();
    CellColumns fs( mesh );
    Field field( fs.createField<int>() );
    EXPECT( field.shape( 0 ) == mesh.cells().size() );
    EXPECT( field.shape( 0 ) == fs.nb_cells() );

    array::ArrayView<int, 1> value     = array::make_view<int, 1>( field );
    array::ArrayView<int, 1> halo      = array::make_view<int, 1>( mesh.cells().halo() );
    array::ArrayView<int, 1> partition = array::make_view<int, 1>( mesh.cells().partition() );
    array::IndexView<idx_t, 1> ridx    = array::make_indexview<idx_t, 1>( mesh.cells().remote_index() );
    array::ArrayView<gidx_t, 1> gidx   = array::make_view<gidx_t, 1>( mesh.cells().global_index() );

    const size_t nb_cells = mesh.cells().size();
    for ( size_t j = 0; j < nb_cells; ++j ) {
        if ( halo( j ) ) {
            value( j ) = -1;
            EXPECT( false );
        }
        else {
            value( j ) = partition( j );
        }
    }

    fs.haloExchange( field );

    for ( size_t j = 0; j < nb_cells; ++j ) {
        EXPECT( value( j ) == partition( j ) );
    }
}

CASE( "test_functionspace_CellColumns_halo_1" ) {
    Mesh mesh = generate_mesh();
    CellColumns fs( mesh, option::halo( 1 ) );

    Field field( fs.createField<double>() );
    EXPECT( field.shape( 0 ) == mesh.cells().size() );
    EXPECT( field.shape( 0 ) == fs.nb_cells() );

    auto value     = array::make_view<double, 1>( field );
    auto halo      = array::make_view<int, 1>( mesh.cells().halo() );
    auto partition = array::make_view<int, 1>( mesh.cells().partition() );
    auto ridx      = array::make_indexview<idx_t, 1>( mesh.cells().remote_index() );
    auto gidx      = array::make_view<gidx_t, 1>( mesh.cells().global_index() );


    const size_t nb_cells = mesh.cells().size();
    for ( size_t j = 0; j < nb_cells; ++j ) {
        if ( halo( j ) ) {
            value( j ) = -1.;
        }
        else {
            value( j ) = partition( j );
        }
    }
    fs.haloExchange( field );
    for ( size_t j = 0; j < nb_cells; ++j ) {
        Log::info() << j << "\t" << gidx( j ) << "\t" << ridx( j ) << "\t" << gidx( ridx( j ) ) << " --- " << value( j )
                    << std::endl;
        EXPECT( value( j ) == partition( j ) );
    }

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
