/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

#include "tests/AtlasTestEnvironment.h"

using namespace std;
using namespace eckit;

//-----------------------------------------------------------------------------

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_rename") {
    FieldSet fieldset;
    auto field_0 = fieldset.add(Field("0", make_datatype<double>(), array::make_shape(10,4)));
    auto field_1 = fieldset.add(Field("1", make_datatype<double>(), array::make_shape(10,5)));
    auto field_2 = fieldset.add(Field("2", make_datatype<double>(), array::make_shape(10,6)));

    EXPECT(fieldset.has("0"));
    EXPECT(fieldset.has("1"));
    EXPECT(fieldset.has("2"));

    field_0.rename("field_0");
    field_1.rename("field_1");
    field_2.rename("field_2");

    EXPECT(fieldset.has("field_0"));
    EXPECT(fieldset.has("field_1"));
    EXPECT(fieldset.has("field_2"));
    EXPECT(!fieldset.has("0"));
    EXPECT(!fieldset.has("1"));
    EXPECT(!fieldset.has("2"));

    EXPECT_EQ(fieldset.field(0).name(),"field_0");
    EXPECT_EQ(fieldset.field(1).name(),"field_1");
    EXPECT_EQ(fieldset.field(2).name(),"field_2");

    EXPECT_EQ(fieldset.field("field_0").shape(1),4);
    EXPECT_EQ(fieldset.field("field_1").shape(1),5);
    EXPECT_EQ(fieldset.field("field_2").shape(1),6);

    field_1.rename("");
    EXPECT(!fieldset.has("field_1"));
    EXPECT_EQ(fieldset.field(1).name(),std::string(""));
    EXPECT(fieldset.has("[1]"));

}

CASE("test_duplicate_name_throws") {
    FieldSet fieldset;
    auto field_0 = fieldset.add(Field("0", make_datatype<double>(), array::make_shape(10,4)));
    auto field_1 = fieldset.add(Field("0", make_datatype<double>(), array::make_shape(10,5))); // same name as field_0, uh-oh !
    auto field_2 = fieldset.add(Field("2", make_datatype<double>(), array::make_shape(10,6)));

    Field f;
    EXPECT_NO_THROW(f = fieldset["2"]); // OK
    EXPECT_EQ(f.shape(1), 6);

    EXPECT_THROWS(f = fieldset["0"]);   // ambiguous because field_0 and field_1 have same name, should throw
    field_1.rename("1");                // fix ambiguity between field_0 and field_1
    EXPECT_NO_THROW(f = fieldset["0"]); // no longer ambiguous
    EXPECT_EQ(f.shape(1), 4);           // to be sure that we got the right field
    EXPECT_NO_THROW(f = fieldset["1"]); // no longer ambiguous
    EXPECT_EQ(f.shape(1), 5);           // to be sure that we got the right field

    field_2.rename("0");                // Introduce new ambiguity between field_0 and field_2
    EXPECT_THROWS(f = fieldset["0"]);   // ambiguous because field_0 and field_2 have same name, should throw
    field_2.rename("2");                // fix ambiguity
    EXPECT_NO_THROW(f = fieldset["0"]); // no longer ambiguous
    EXPECT_EQ(f.shape(1), 4);           // to be sure we got the right field
    EXPECT_NO_THROW(f = fieldset["2"]); // no longer ambiguous
    EXPECT_EQ(f.shape(1), 6);           // to be sure we got the right field
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
