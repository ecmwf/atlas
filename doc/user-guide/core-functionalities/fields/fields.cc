#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Metadata.h"

using atlas::Field;
using atlas::FieldSet;
using atlas::Log;
using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::array::make_view;

int main( int argc, char* argv[] ) {
    atlas::initialise( argc, argv );

    // Define fields
    Field field_pressure( "pressure", make_datatype<double>(), make_shape( 100 ) );
    Field field_wind( "wind", make_datatype<double>(), make_shape( 100, 2 ) );

    // Access field data
    auto pressure = make_view<double, 1>( field_pressure );
    auto wind     = make_view<double, 2>( field_wind );

    // Assign values to fields
    for ( size_t jnode = 0; jnode < 100; ++jnode ) {
        pressure( jnode )          = 101325.0;
        wind( jnode, size_t( 0 ) ) = 0.01 + double( jnode );
        wind( jnode, size_t( 1 ) ) = 0.02 + double( jnode );
    }

    // Add info to fields
    std::string unitsP, unitsW;
    field_pressure.metadata().set( "units", "[Pa]" );
    field_pressure.metadata().get( "units", unitsP );
    field_wind.metadata().set( "units", "[m/s]" );
    field_wind.metadata().get( "units", unitsW );

    // Define fieldSet
    FieldSet fields;
    fields.add( field_pressure );  // Add field_pressure to fieldSet
    fields.add( field_wind );      // Add field_wind to fieldSet

    // Retrieve field from fieldSet
    Field field_pressure2 = fields.field( "pressure" );
    Field field_wind2     = fields.field( "wind" );

    // Print some useful info
    Log::info() << "name   = " << field_wind.name() << std::endl;
    Log::info() << "size   = " << field_wind.size() << std::endl;
    Log::info() << "units  = " << unitsW << std::endl;
    Log::info() << "rank   = " << field_wind.rank() << std::endl;
    Log::info() << "shape  = " << field_wind.shape( 0 ) << "    " << field_wind.shape( 1 ) << std::endl;
    Log::info() << "memory = " << field_wind.bytes() << " bytes" << std::endl;
    Log::info() << "type   = " << field_wind.datatype().str() << std::endl;
    Log::info() << "kind   = " << field_wind.datatype().kind() << std::endl;

    // Print some values
    Log::info() << "pressure(9) = " << pressure( 9 ) << std::endl;
    Log::info() << "wind(9, 0)  = " << wind( 9, 0 ) << std::endl;
    Log::info() << "wind(9, 1)  = " << wind( 9, 1 ) << std::endl;

    atlas::finalize();
    atlas::mpi::finalize();

    return 0;
}
