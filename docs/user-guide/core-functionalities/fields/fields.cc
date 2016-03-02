#include "atlas/atlas.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/util/Metadata.h"

using namespace std;
using namespace atlas;
using namespace atlas::field;
using namespace atlas::array;

int main(int argc, char *argv[])
{
    atlas_init(argc, argv);

    // Define fields
    Field::Ptr pressureField(Field::create<double>
                             ("pressure", make_shape(100)));
    Field::Ptr windField(Field::create<double>
                        ("wind", make_shape(2, 100)));

    // Initialize fields
    ArrayView <double,1> pressure(*pressureField);
    ArrayView <double,2> wind    (*windField);

    // Assign values to fields
    for (int jnode = 0; jnode < 100; ++jnode)
    {
        pressure(jnode) = 101325.0;
        wind(0,jnode)   = 0.01 + double(jnode);
        wind(1,jnode)   = 0.02 + double(jnode);
    }

    // Add info to fields
    string unitsP, unitsW;
    pressureField->metadata().set("units", "[Pa]");
    pressureField->metadata().get("units", unitsP);
    windField    ->metadata().set("units", "[m/s]");
    windField    ->metadata().get("units", unitsW);

    // Define fieldSet
    FieldSet fields;
    fields.add(*pressureField); // Add pressureField to fieldSet
    fields.add(*windField);     // Add windField to fieldSet

    // Retrieve field from fieldSet
    Field& pressureField2 = fields.field("pressure");
    Field& windField2     = fields.field("wind");

    // Print some useful info
    cout << "name   = " << windField->name()     << endl;
    cout << "size   = " << windField->size()     << endl;
    cout << "units  = " << unitsW                << endl;
    cout << "rank   = " << windField->rank()     << endl;
    cout << "shape  = " << windField->shape(1)   << "    "
                        << windField->shape(2)   << endl;
    cout << "memory = " << windField->bytes()
                        << " bytes"              << endl;
    cout << "type   = " << windField->datatype().str()  << endl;
    cout << "kind   = " << windField->datatype().kind() << endl;

    // Print some values
    cout << "pressure(9) = " << pressure(9) << endl;
    cout << "wind(0, 9)  = " << wind(0,9)   << endl;
    cout << "wind(1, 9)  = " << wind(1,9)   << endl;

    atlas_finalize();

    return 0;
}
