/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <string>

#include "eckit/log/JSON.h"
#include "eckit/parser/JSONParser.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/DataType.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/MultiField.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::field;
using namespace atlas::field;

namespace atlas {
namespace test {

// -------------------------------------------------------------------
// Example IFS MultiField creato

// ---  Declaration (in .h file)
class MultiFieldCreatorIFS : public MultiFieldCreator {
public:
    MultiFieldCreatorIFS(const eckit::Parametrisation& p = util::Config()): MultiFieldCreator(p) {}
    ~MultiFieldCreatorIFS() override = default;
    void generate(MultiField& fieldpool, const eckit::Parametrisation& p = util::Config()) const override;
};

// ---  Implementation (in .cc file)
void MultiFieldCreatorIFS::generate(MultiField& multifield, const eckit::Parametrisation& p) const {
    const eckit::LocalConfiguration* params = dynamic_cast<const eckit::LocalConfiguration*>(&p);
    if (!params) {
        throw_Exception("Parametrisation has to be of atlas::util::Config type");
    }


    long nproma = 0;
    long ngptot = 0;
    ;
    long nfld = 0;
    ;
    long nlev = 0;
    ;
    long nblk = 0;
    ;

    if (!p.get("ngptot", ngptot)) {
        std::string grid_uid;
        if (p.get("grid", grid_uid)) {
            ngptot = Grid(grid_uid).size();
        }
        else {
            throw_Exception("Could not find 'ngptot' in Parametrisation");
        }
    }

    if (!p.get("nproma", nproma)) {
        throw_Exception("Could not find 'nproma' in Parametrisation");
    }

    if (!p.get("nlev", nlev)) {
        throw_Exception("Could not find 'nlev' in Parametrisation");
    }

    array::DataType datatype = array::DataType::create<double>();
    std::string datatype_str;
    if (p.get("datatype", datatype_str)) {
        datatype = array::DataType(datatype_str);
    }
    else {
        array::DataType::kind_t kind(array::DataType::kind<double>());
        p.get("kind", kind);
        if (!array::DataType::kind_valid(kind)) {
            std::stringstream msg;
            msg << "Could not create field. kind parameter unrecognized";
            throw_Exception(msg.str());
        }
        datatype = array::DataType(kind);
    }

    nblk = std::ceil(static_cast<double>(ngptot) / static_cast<double>(nproma));

    std::vector<eckit::LocalConfiguration> fields;
    params->get("fields", fields);
    nfld = fields.size();

    auto multiarray_spec = array::make_shape(nblk, nfld, nlev, nproma);
    auto& multiarray     = multifield.allocate(datatype, multiarray_spec);

    auto field_shape      = array::make_shape(multiarray.shape(0), multiarray.shape(2), multiarray.shape(3));
    auto field_strides    = array::make_strides(multiarray.stride(0), multiarray.stride(2), multiarray.stride(3));
    auto field_array_spec = array::ArraySpec(field_shape, field_strides);

    for (size_t i = 0; i < fields.size(); ++i) {
        std::string name;
        fields[i].get("name", name);
        Field field;
        constexpr auto all = array::Range::all();
        if (datatype.kind() == array::DataType::KIND_REAL64) {
            auto slice = array::make_view<double, 4>(multiarray).slice(all, i, all, all);
            field      = Field(name, slice.data(), field_array_spec);
        }
        else if (datatype.kind() == array::DataType::KIND_REAL32) {
            auto slice = array::make_view<float, 4>(multiarray).slice(all, i, all, all);
            field      = Field(name, slice.data(), field_array_spec);
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }
        field.set_levels(nlev);
        multifield.fields_.emplace_back(field);
        multifield.field_index_[field.name()] = i;
    }
}

// Register in factory
MultiFieldCreatorBuilder<MultiFieldCreatorIFS> __MultiFieldCreatorIFS("MultiFieldCreatorIFS");

// ===================================================================
//                               BEGIN TESTS
// ===================================================================


CASE("multifield_generator") {
    EXPECT(MultiFieldCreatorFactory::has("MultiFieldCreatorIFS"));
    std::unique_ptr<MultiFieldCreator> MultiFieldCreator(MultiFieldCreatorFactory::build("MultiFieldCreatorIFS"));
}


CASE("multifield_create") {
    int nproma = 16;
    int nlev   = 100;
    int ngptot = 2000;

    auto json = [&]() -> std::string {
        util::Config p;
        p.set("ngptot", ngptot);
        p.set("nproma", nproma);
        p.set("nlev", nlev);
        p.set("datatype", array::DataType::real64().str());

        std::vector<util::Config> fields(3);
        fields[0].set("name", "temperature");
        fields[1].set("name", "pressure");
        fields[2].set("name", "density");

        p.set("fields", fields);

        // We can also translate parameters to a json:
        std::stringstream json;
        eckit::JSON js(json);
        js << p;
        std::string json_str = json.str();
        Log::info() << "json = " << json_str << std::endl;
        return json_str;
    };


    // And we can create back parameters from json:
    std::stringstream json_stream;
    json_stream << json();
    util::Config config(json_stream);

    MultiField multifield("MultiFieldCreatorIFS", config);

    EXPECT(multifield.has("temperature"));
    EXPECT(multifield.has("pressure"));
    EXPECT(multifield.has("density"));

    Log::info() << multifield.field("temperature") << std::endl;
    Log::info() << multifield.field("pressure") << std::endl;
    Log::info() << multifield.field("density") << std::endl;

    auto temp = array::make_view<double, 3>(multifield.field("temperature"));
    auto pres = array::make_view<double, 3>(multifield.field("pressure"));
    auto dens = array::make_view<double, 3>(multifield.field("density"));

    temp(1, 2, 3)   = 4;
    pres(5, 6, 7)   = 8;
    dens(9, 10, 11) = 12;

    EXPECT_EQ(temp.stride(2), 1);
    EXPECT_EQ(temp.stride(1), nproma);
    EXPECT_EQ(temp.stride(0), nproma * nlev * multifield.size());

    // Underlying array of MultiField
    {
        auto multiarray = array::make_view<double, 4>(multifield);
        EXPECT_EQ(multiarray.stride(3), 1);
        EXPECT_EQ(multiarray.stride(2), nproma);
        EXPECT_EQ(multiarray.stride(1), nproma * nlev);
        EXPECT_EQ(multiarray.stride(0), nproma * nlev * multifield.size());

        EXPECT_EQ(multiarray(1, 0, 2, 3), 4.);
        EXPECT_EQ(multiarray(5, 1, 6, 7), 8.);
        EXPECT_EQ(multiarray(9, 2, 10, 11), 12.);
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
