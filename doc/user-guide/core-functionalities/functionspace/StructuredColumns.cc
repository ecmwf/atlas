#include <cmath>
#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"

using namespace atlas;
using atlas::StructuredMeshGenerator;
using atlas::array::make_view;
using atlas::functionspace::StructuredColumns;
using atlas::output::Gmsh;

int main( int argc, char* argv[] ) {
    atlas::initialise( argc, argv );

    // Generate global reduced grid
    Grid grid( "N32" );

    // Generate functionspace associated to grid
    StructuredColumns fs( grid );

    // Variables for scalar1 field definition
    const double deg2rad = M_PI / 180.;
    const double zlatc   = 0.0 * M_PI;
    const double zlonc   = 1.0 * M_PI;
    const double zrad    = 2.0 * M_PI / 9.0;
    int jnode            = 0;

    // Calculate scalar function
    Field field_scalar1 = fs.createField<double>( option::name( "scalar1" ) );
    auto xy             = make_view<double, 2>( fs.xy() );
    auto scalar1        = make_view<double, 1>( field_scalar1 );

    for ( idx_t j = fs.j_begin(); j < fs.j_end(); ++j ) {
        for ( idx_t i = fs.i_begin( j ); i < fs.i_end( j ); ++i ) {
            double zlon  = xy( fs.index( i, j ), XX ) * deg2rad;
            double zlat  = xy( fs.index( i, j ), YY ) * deg2rad;
            double zdist = 2.0 * std::sqrt( ( cos( zlat ) * std::sin( ( zlon - zlonc ) / 2. ) ) *
                                                ( std::cos( zlat ) * std::sin( ( zlon - zlonc ) / 2. ) ) +
                                            std::sin( ( zlat - zlatc ) / 2. ) * std::sin( ( zlat - zlatc ) / 2. ) );

            scalar1( jnode ) = 0.0;
            if ( zdist < zrad ) {
                scalar1( jnode ) = 0.5 * ( 1. + std::cos( M_PI * zdist / zrad ) );
            }
            ++jnode;
        }
    }

    // Write field
    {
        // Generate visualisation mesh associated to grid
        StructuredMeshGenerator meshgenerator;
        Mesh mesh = meshgenerator.generate( grid );

        Gmsh gmsh( "scalar1.msh" );
        gmsh.write( mesh );
        gmsh.write( field_scalar1 );
    }

    atlas::finalise();
    atlas::mpi::finalise();
    return 0;
}
