#include "atlas/library/atlas.h"
#include "atlas/runtime/Log.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/util/Metadata.h"

using atlas::atlas_init;
using atlas::atlas_finalize;
using atlas::Log;
using atlas::field::Field;
using atlas::field::FieldSet;
using atlas::array::ArrayView;
using atlas::array::make_shape;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    // Define fields
    Field::Ptr field_pressure(
      Field::create<double>("pressure", make_shape(100)) );
    Field::Ptr field_wind(
      Field::create<double>("wind", make_shape(100, 2)) );

    // Access field data
    ArrayView <double,1> pressure(*field_pressure);
    ArrayView <double,2> wind    (*field_wind);

    // Assign values to fields
    for (size_t jnode = 0; jnode < 100; ++jnode)
    {
        pressure(jnode) = 101325.0;
        wind(jnode,0)   = 0.01 + double(jnode);
        wind(jnode,1)   = 0.02 + double(jnode);
    }

    // Add info to fields
    std::string unitsP, unitsW;
    field_pressure->metadata().set("units", "[Pa]");
    field_pressure->metadata().get("units", unitsP);
    field_wind    ->metadata().set("units", "[m/s]");
    field_wind    ->metadata().get("units", unitsW);

    // Define fieldSet
    FieldSet fields;
    fields.add(*field_pressure); // Add field_pressure to fieldSet
    fields.add(*field_wind);     // Add field_wind to fieldSet

    // Retrieve field from fieldSet
    Field& field_pressure2 = fields.field("pressure");
    Field& field_wind2     = fields.field("wind");

    // Print some useful info
    Log::info() << "name   = " << field_wind->name()     << std::endl;
    Log::info() << "size   = " << field_wind->size()     << std::endl;
    Log::info() << "units  = " << unitsW                 << std::endl;
    Log::info() << "rank   = " << field_wind->rank()     << std::endl;
    Log::info() << "shape  = " << field_wind->shape(0)   << "    "
                               << field_wind->shape(1)   << std::endl;
    Log::info() << "memory = " << field_wind->bytes()
                               << " bytes"              << std::endl;
    Log::info() << "type   = " << field_wind->datatype().str()  << std::endl;
    Log::info() << "kind   = " << field_wind->datatype().kind() << std::endl;

    // Print some values
    Log::info() << "pressure(9) = " << pressure(9) << std::endl;
    Log::info() << "wind(9, 0)  = " << wind(9,0)   << std::endl;
    Log::info() << "wind(9, 1)  = " << wind(9,1)   << std::endl;

    atlas_finalize();

    return 0;
}
