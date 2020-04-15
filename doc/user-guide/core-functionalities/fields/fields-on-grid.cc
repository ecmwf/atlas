#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

using atlas::Field;
using atlas::Log;
using atlas::StructuredGrid;
using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::array::make_view;

int main( int argc, char* argv[] ) {
    atlas::initialise( argc, argv );

    int jnode            = 0;
    const double rpi     = 2.0 * asin( 1.0 );
    const double deg2rad = rpi / 180.;
    const double zlatc   = 0.0 * rpi;
    const double zlonc   = 1.0 * rpi;
    const double zrad    = 2.0 * rpi / 9.0;
    double zdist, zlon, zlat;

    StructuredGrid grid( "N32" );
    const size_t nb_nodes = grid.size();

    Field field_pressure( "pressure", make_datatype<double>(), make_shape( nb_nodes ) );

    auto pressure = make_view<double, 1>( field_pressure );
    for ( size_t jlat = 0; jlat < grid.ny(); ++jlat ) {
        zlat = grid.y( jlat );
        zlat = zlat * deg2rad;
        for ( size_t jlon = 0; jlon < grid.nx( jlat ); ++jlon ) {
            zlon  = grid.x( jlon, jlat );
            zlon  = zlon * deg2rad;
            zdist = 2.0 * sqrt( ( cos( zlat ) * sin( ( zlon - zlonc ) / 2 ) ) *
                                    ( cos( zlat ) * sin( ( zlon - zlonc ) / 2 ) ) +
                                sin( ( zlat - zlatc ) / 2 ) * sin( ( zlat - zlatc ) / 2 ) );

            pressure( jnode ) = 0.0;
            if ( zdist < zrad ) {
                pressure( jnode ) = 0.5 * ( 1. + cos( rpi * zdist / zrad ) );
            }
            jnode = jnode + 1;
        }
    }

    Log::info() << "==========================================" << std::endl;
    Log::info() << "memory field_pressure = " << field_pressure.bytes() * 1.e-9 << " GB" << std::endl;
    Log::info() << "==========================================" << std::endl;

    atlas::finalise();
    atlas::mpi::finalize();

    return 0;
}
