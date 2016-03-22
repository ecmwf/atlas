#include "atlas/atlas.h"
#include "atlas/runtime/Log.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/util/Metadata.h"

using namespace std;
using namespace eckit;
using namespace atlas;
using namespace atlas::field;
using namespace atlas::array;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    // Define fields
    SharedPtr<Field> field_pressure(
      Field::create<double>("pressure", make_shape(100)) );
    SharedPtr<Field> field_wind(
      Field::create<double>("wind", make_shape(100, 2)) );

    // Access field data
    ArrayView <double,1> pressure(*field_pressure);
    ArrayView <double,2> wind    (*field_wind);

    // Assign values to fields
    for (int jnode = 0; jnode < 100; ++jnode)
    {
        pressure(jnode) = 101325.0;
        wind(jnode,0)   = 0.01 + double(jnode);
        wind(jnode,1)   = 0.02 + double(jnode);
    }

    // Add info to fields
    string unitsP, unitsW;
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
    Log::info() << "name   = " << field_wind->name()     << endl;
    Log::info() << "size   = " << field_wind->size()     << endl;
    Log::info() << "units  = " << unitsW                << endl;
    Log::info() << "rank   = " << field_wind->rank()     << endl;
    Log::info() << "shape  = " << field_wind->shape(0)   << "    "
                               << field_wind->shape(1)   << endl;
    Log::info() << "memory = " << field_wind->bytes()
                               << " bytes"              << endl;
    Log::info() << "type   = " << field_wind->datatype().str()  << endl;
    Log::info() << "kind   = " << field_wind->datatype().kind() << endl;

    // Print some values
    Log::info() << "pressure(9) = " << pressure(9) << endl;
    Log::info() << "wind(9, 0)  = " << wind(9,0)   << endl;
    Log::info() << "wind(9, 1)  = " << wind(9,1)   << endl;

    atlas_finalize();

    return 0;
}
