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

CASE("test_fieldset_concatenation") {
    FieldSet fieldset_1;
    FieldSet fieldset_2;
    auto field_1 = fieldset_1.add(Field("0", make_datatype<double>(), array::make_shape(10,4)));
    auto field_1_v = array::make_view<double,2>(field_1);
    field_1_v(1,1) = 2.;

    fieldset_2.add_fieldset(fieldset_1);
    auto field_2 = fieldset_2.field("0");
    field_2.rename("");
    auto field_2_v = array::make_view<double,2>(field_2);
    field_2_v(1,1) = 1.;

    EXPECT(!fieldset_2.has(""));
    EXPECT_EQ(fieldset_2.field(0).name(), "");
    EXPECT_EQ(field_1_v(1,1), 1.);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
