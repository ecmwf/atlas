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

#include "eckit/config/YAMLConfiguration.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/DataType.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/MultiField.h"
#include "atlas/grid/Grid.h"
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
    MultiFieldCreatorIFS(const eckit::Configuration& config = util::Config()): MultiFieldCreator(config) {}
    ~MultiFieldCreatorIFS() override = default;
    MultiFieldImpl* create(const eckit::Configuration& = util::Config()) const override;
};

// ---  Implementation (in .cc file)
MultiFieldImpl* MultiFieldCreatorIFS::create(const eckit::Configuration& config) const {
    long ngptot = config.getLong("ngptot");
    long nproma = config.getLong("nproma");
    long nlev   = config.getLong("nlev");
    long nblk   = 0;


    array::DataType datatype = array::DataType::create<double>();
    std::string datatype_str;
    if (config.get("datatype", datatype_str)) {
        datatype = array::DataType(datatype_str);
    }
    else {
        array::DataType::kind_t kind(array::DataType::kind<double>());
        config.get("kind", kind);
        if (!array::DataType::kind_valid(kind)) {
            std::stringstream msg;
            msg << "Could not create field. kind parameter unrecognized";
            throw_Exception(msg.str());
        }
        datatype = array::DataType(kind);
    }

    nblk = std::ceil(static_cast<double>(ngptot) / static_cast<double>(nproma));

    auto fields = config.getSubConfigurations("fields");
    long nfld   = 0;
    for (const auto& field_config : fields) {
        long nvar = 1;
        field_config.get("nvar", nvar);
        nfld += nvar;
    }

    auto multiarray_shape = array::make_shape(nblk, nfld, nlev, nproma);

    MultiFieldImpl* multifield = new MultiFieldImpl{array::ArraySpec{datatype, multiarray_shape}};

    auto& multiarray = multifield->array();

    size_t multiarray_field_idx = 0;
    for (size_t i = 0; i < fields.size(); ++i) {
        std::string name;
        fields[i].get("name", name);
        Field field;
        size_t field_vars = 1;

        if (fields[i].get("nvar", field_vars)) {
            auto field_shape =
                array::make_shape(multiarray.shape(0), field_vars, multiarray.shape(2), multiarray.shape(3));
            auto field_strides    = multiarray.strides();
            auto field_array_spec = array::ArraySpec(field_shape, field_strides);

            constexpr auto all = array::Range::all();
            const auto range   = array::Range(multiarray_field_idx, multiarray_field_idx + field_vars);
            if (datatype.kind() == array::DataType::KIND_REAL64) {
                auto slice = array::make_view<double, 4>(multiarray).slice(all, range, all, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else if (datatype.kind() == array::DataType::KIND_REAL32) {
                auto slice = array::make_view<float, 4>(multiarray).slice(all, range, all, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
            field.set_variables(field_vars);
        }
        else {
            auto field_shape   = array::make_shape(multiarray.shape(0), multiarray.shape(2), multiarray.shape(3));
            auto field_strides = array::make_strides(multiarray.stride(0), multiarray.stride(2), multiarray.stride(3));
            auto field_array_spec = array::ArraySpec(field_shape, field_strides);

            constexpr auto all = array::Range::all();
            if (datatype.kind() == array::DataType::KIND_REAL64) {
                auto slice = array::make_view<double, 4>(multiarray).slice(all, multiarray_field_idx, all, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else if (datatype.kind() == array::DataType::KIND_REAL32) {
                auto slice = array::make_view<float, 4>(multiarray).slice(all, multiarray_field_idx, all, all);
                field      = Field(name, slice.data(), field_array_spec);
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
        }
        field.set_levels(nlev);

        multifield->add(field);

        multiarray_field_idx += field_vars;
    }
    return multifield;
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
    using Value = float;
    int nproma  = 16;
    int nlev    = 100;
    int ngptot  = 2000;

    auto json = [&]() -> std::string {
        util::Config p;
        p.set("type", "MultiFieldCreatorIFS");
        p.set("ngptot", ngptot);
        p.set("nproma", nproma);
        p.set("nlev", nlev);
        p.set("datatype", array::make_datatype<Value>().str());
        p.set("fields", {
                            util::Config("name", "temperature"),     //
                            util::Config("name", "pressure"),        //
                            util::Config("name", "density"),         //
                            util::Config("name", "clv")("nvar", 5),  //
                            util::Config("name", "wind_u")           //
                        });
        return p.json();
    };

    SECTION("Print configuration") {
        Log::info() << "json = " << json() << std::endl;
    }

    SECTION("test") {
        MultiField multifield{eckit::YAMLConfiguration{json()}};

        const auto nblk = multifield.array().shape(0);
        const auto nvar = multifield.array().shape(1);
        const auto nfld = multifield.size();
        EXPECT_EQ(nfld, 5);
        EXPECT_EQ(nvar, 9);

        EXPECT_EQ(multifield.size(), 5);
        EXPECT(multifield.has("temperature"));
        EXPECT(multifield.has("pressure"));
        EXPECT(multifield.has("density"));
        EXPECT(multifield.has("clv"));
        EXPECT(multifield.has("wind_u"));

        Log::info() << multifield.field("temperature") << std::endl;
        Log::info() << multifield.field("pressure") << std::endl;
        Log::info() << multifield.field("density") << std::endl;
        Log::info() << multifield.field("clv") << std::endl;
        Log::info() << multifield.field("wind_u") << std::endl;

        auto temp   = array::make_view<Value, 3>(multifield.field("temperature"));
        auto pres   = array::make_view<Value, 3>(multifield.field("pressure"));
        auto dens   = array::make_view<Value, 3>(multifield.field("density"));
        auto clv    = array::make_view<Value, 4>(multifield.field("clv"));  // note rank 4
        auto wind_u = array::make_view<Value, 3>(multifield.field("wind_u"));

        EXPECT_EQ(multifield[0].name(), "temperature");
        EXPECT_EQ(multifield[1].name(), "pressure");
        EXPECT_EQ(multifield[2].name(), "density");
        EXPECT_EQ(multifield[3].name(), "clv");
        EXPECT_EQ(multifield[4].name(), "wind_u");

        auto block_stride  = multifield.array().stride(0);
        auto field_stride  = nproma * nlev;
        auto level_stride  = nproma;
        auto nproma_stride = 1;

        temp(1, 2, 3)      = 4;
        pres(5, 6, 7)      = 8;
        dens(9, 10, 11)    = 12;
        clv(13, 2, 14, 15) = 16;
        wind_u(17, 18, 3)  = 19;

        EXPECT_EQ(temp.stride(0), block_stride);
        EXPECT_EQ(temp.stride(1), level_stride);
        EXPECT_EQ(temp.stride(2), nproma_stride);
        EXPECT_EQ(temp.size(), nblk * nlev * nproma);

        EXPECT_EQ(clv.stride(0), block_stride);
        EXPECT_EQ(clv.stride(1), field_stride);
        EXPECT_EQ(clv.stride(2), level_stride);
        EXPECT_EQ(clv.stride(3), nproma_stride);

        EXPECT_EQ(clv.size(), nblk * 5 * nlev * nproma);


        // Advanced usage, to access underlying array. This should only be used
        // in a driver and not be exposed to algorithms.
        {
            auto multiarray = array::make_view<Value, 4>(multifield);
            EXPECT_EQ(multiarray.stride(0), block_stride);
            EXPECT_EQ(multiarray.stride(1), field_stride);
            EXPECT_EQ(multiarray.stride(2), level_stride);
            EXPECT_EQ(multiarray.stride(3), nproma_stride);

            EXPECT_EQ(multiarray(1, 0, 2, 3), 4.);
            EXPECT_EQ(multiarray(5, 1, 6, 7), 8.);
            EXPECT_EQ(multiarray(9, 2, 10, 11), 12.);
            EXPECT_EQ(multiarray(13, 5, 14, 15), 16.);
            EXPECT_EQ(multiarray(17, 8, 18, 3), 19.);

            EXPECT_EQ(multiarray.size(), nblk * nvar * nlev * nproma);
        }

        // access FieldSet through MultiField
        auto fieldset = multifield->fieldset();
        auto field_v = array::make_view<Value,3>(fieldset.field("temperature"));
        EXPECT_EQ(fieldset.size(), 5);
        EXPECT(fieldset.has("temperature"));
        EXPECT(fieldset.has("wind_u"));
        EXPECT_EQ(field_v(1,2,3), 4);
    }

    SECTION("test registry") {
        {
            Field field = MultiField {eckit::YAMLConfiguration{json()}}.field("temperature");
            auto temp = array::make_view<Value,3>(field);
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
