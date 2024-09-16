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
#include <vector>

#include "eckit/config/YAMLConfiguration.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/MultiField.h"
#include "atlas/field/MultiFieldCreatorArray.h"
#include "atlas/field/detail/MultiFieldImpl.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::field;

namespace atlas {
namespace test {

const std::vector<int> make_shape(const std::initializer_list<int>& list) {
    return std::vector<int>(list);
}

CASE("multifield_generator") {
    EXPECT(MultiFieldCreatorFactory::has("MultiFieldCreatorIFS"));
    std::unique_ptr<MultiFieldCreator> MultiFieldCreatorIFS(MultiFieldCreatorFactory::build("MultiFieldCreatorIFS"));

    EXPECT(MultiFieldCreatorFactory::has("MultiFieldCreatorArray"));
    std::unique_ptr<MultiFieldCreator> MultiFieldCreatorArray(MultiFieldCreatorFactory::build("MultiFieldCreatorArray"));
}


CASE("multifield_ifs_create") {
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
                            util::Config("name", "temperature"),
                            util::Config("name", "pressure"),
                            util::Config("name", "density"),
                            util::Config("name", "clv")("nvar", 5),  // 'clv' with 5 subvariables
                            util::Config("name", "wind_u")
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

CASE("multifield_array_create") {
    using Value = float;
    int nproma  = 16;
    int nlev    = 100;
    int ngptot  = 2000;

    const int nblks = (ngptot + nproma - 1) / nproma;
    const std::vector<std::string> var_names = {"temperature", "pressure", "density", "clv", "wind_u"};

    auto json = [&]() -> std::string {
        util::Config p;
        p.set("type", "MultiFieldCreatorArray");
        p.set("shape", {nblks, -1, nlev, nproma});
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

    SECTION("test_MultiFieldArray_noconfig_3d") {
        int nlev = 3;
        const std::vector<int> vshape = make_shape({nblks, -1, nlev, nproma});
        MultiField multifield(array::make_datatype<Value>(), vshape, var_names);

        const auto nblk = multifield.array().shape(0);
        const auto nvar = multifield.array().shape(1);
        nlev = multifield.array().shape(2);
        const auto nfld = multifield.size();
        EXPECT_EQ(nfld, 5);
        EXPECT_EQ(nvar, 5);

        EXPECT_EQ(multifield.size(), 5);
        EXPECT(multifield.has("temperature"));
        EXPECT(multifield.has("clv"));

        Log::info() << multifield.field("temperature") << std::endl;
        Log::info() << multifield.field("clv") << std::endl;

        auto temp   = array::make_view<Value, 3>(multifield.field("temperature"));
        auto clv    = array::make_view<Value, 3>(multifield.field("clv"));
        EXPECT_EQ(multifield[0].name(), "temperature");
        EXPECT_EQ(multifield[3].name(), "clv");

        auto block_stride  = multifield.array().stride(0);
        auto field_stride  = nproma * nlev;
        auto level_stride  = nproma;
        auto nproma_stride = 1;

        temp(1, 2, 3)      = 4;
        clv(13, 2, 14)     = 16;

        EXPECT_EQ(temp.stride(0), block_stride);
        EXPECT_EQ(temp.stride(1), level_stride);
        EXPECT_EQ(temp.stride(2), nproma_stride);
        EXPECT_EQ(temp.size(), nblk * nlev * nproma);

        // Advanced usage, to access underlying array. This should only be used
        // in a driver and not be exposed to algorithms.
        {
            auto multiarray = array::make_view<Value, 4>(multifield);
            EXPECT_EQ(multiarray.stride(0), block_stride);
            EXPECT_EQ(multiarray.stride(1), field_stride);
            EXPECT_EQ(multiarray.stride(2), level_stride);
            EXPECT_EQ(multiarray.stride(3), nproma_stride);

            EXPECT_EQ(multiarray(1, 0, 2, 3), 4.);
            EXPECT_EQ(multiarray(13, 3, 2, 14), 16.);

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

    SECTION("test_MultiFieldArray_noconfig_2d") {
        const std::vector<int> vshape = make_shape({nblks, -1, nproma});
        MultiField multifield(array::make_datatype<Value>(), vshape, var_names);

        const auto nblk = multifield.array().shape(0);
        const auto nvar = multifield.array().shape(1);
        nlev = multifield.array().shape(2);
        const auto nfld = multifield.size();
        EXPECT_EQ(nfld, 5);
        EXPECT_EQ(nvar, 5);

        EXPECT_EQ(multifield.size(), 5);
        EXPECT(multifield.has("temperature"));
        EXPECT(multifield.has("clv"));

        Log::info() << multifield.field("temperature") << std::endl;
        Log::info() << multifield.field("clv") << std::endl;

        auto temp   = array::make_view<Value, 2>(multifield.field("temperature"));
        auto clv    = array::make_view<Value, 2>(multifield.field("clv"));
        EXPECT_EQ(multifield[0].name(), "temperature");
        EXPECT_EQ(multifield[3].name(), "clv");

        auto block_stride  = multifield.array().stride(0);
        auto field_stride  = nproma;
        auto nproma_stride = 1;

        temp(1, 3)      = 4;
        clv(13, 14)     = 16;

        EXPECT_EQ(temp.stride(0), block_stride);
        EXPECT_EQ(temp.stride(1), nproma_stride);
        EXPECT_EQ(temp.size(), nblk * nproma);

        // Advanced usage, to access underlying array. This should only be used
        // in a driver and not be exposed to algorithms.
        {
            auto multiarray = array::make_view<Value, 3>(multifield);
            EXPECT_EQ(multiarray.stride(0), block_stride);
            EXPECT_EQ(multiarray.stride(1), field_stride);
            EXPECT_EQ(multiarray.stride(2), nproma_stride);

            EXPECT_EQ(multiarray(1, 0, 3), 4.);
            EXPECT_EQ(multiarray(13, 3, 14), 16.);

            EXPECT_EQ(multiarray.size(), nblk * nvar * nproma);
        }

        // access FieldSet through MultiField
        auto fieldset = multifield->fieldset();
        auto field_v = array::make_view<Value, 2>(fieldset.field("temperature"));
        EXPECT_EQ(fieldset.size(), 5);
        EXPECT(fieldset.has("temperature"));
        EXPECT(fieldset.has("wind_u"));
        EXPECT_EQ(field_v(1,3), 4);
    }

    SECTION("Print configuration") {
        Log::info() << "json = " << json() << std::endl;
    }

    SECTION("test_MultiFieldArray_config") {
        MultiField multifield{eckit::YAMLConfiguration{json()}};

        const auto nblk = multifield.array().shape(0);
        const auto nvar = multifield.array().shape(1);
        const auto nfld = multifield.size();
        EXPECT_EQ(nfld, 9);
        EXPECT_EQ(nvar, 9);

        EXPECT_EQ(multifield.size(), 9);
        EXPECT(multifield.has("temperature"));
        EXPECT(multifield.has("pressure"));
        EXPECT(multifield.has("density"));
        EXPECT(multifield.has("clv_0"));
        EXPECT(multifield.has("wind_u"));

        Log::info() << multifield.field("temperature") << std::endl;
        Log::info() << multifield.field("pressure") << std::endl;
        Log::info() << multifield.field("density") << std::endl;
        Log::info() << multifield.field("clv_0") << std::endl;
        Log::info() << multifield.field("wind_u") << std::endl;

        auto temp   = array::make_view<Value, 3>(multifield.field("temperature"));
        auto pres   = array::make_view<Value, 3>(multifield.field("pressure"));
        auto dens   = array::make_view<Value, 3>(multifield.field("density"));
        auto clv    = array::make_view<Value, 3>(multifield.field("clv_0"));  // note rank 4
        auto wind_u = array::make_view<Value, 3>(multifield.field("wind_u"));

        EXPECT_EQ(multifield[0].name(), "temperature");
        EXPECT_EQ(multifield[1].name(), "pressure");
        EXPECT_EQ(multifield[2].name(), "density");
        EXPECT_EQ(multifield[3].name(), "clv_0");
        EXPECT_EQ(multifield[8].name(), "wind_u");

        auto block_stride  = multifield.array().stride(0);
        auto field_stride  = nproma * nlev;
        auto level_stride  = nproma;
        auto nproma_stride = 1;

        temp(1, 2, 3)      = 4;
        pres(5, 6, 7)      = 8;
        dens(9, 10, 11)    = 12;
        clv(13, 14, 15) = 16;
        wind_u(17, 18, 3)  = 19;

        EXPECT_EQ(temp.stride(0), block_stride);
        EXPECT_EQ(temp.stride(1), level_stride);
        EXPECT_EQ(temp.stride(2), nproma_stride);
        EXPECT_EQ(temp.size(), nblk * nlev * nproma);

        EXPECT_EQ(clv.stride(0), block_stride);
        EXPECT_EQ(clv.stride(1), level_stride);
        EXPECT_EQ(clv.stride(2), nproma_stride);

        EXPECT_EQ(clv.size(), nblk * nlev * nproma);


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
            EXPECT_EQ(multiarray(13, 3, 14, 15), 16.);
            EXPECT_EQ(multiarray(17, 8, 18, 3), 19.);

            EXPECT_EQ(multiarray.size(), nblk * nvar * nlev * nproma);
        }

        // access FieldSet through MultiField
        auto fieldset = multifield->fieldset();
        auto field_v = array::make_view<Value,3>(fieldset.field("temperature"));
        EXPECT_EQ(fieldset.size(), 9);
        EXPECT(fieldset.has("temperature"));
        EXPECT(fieldset.has("wind_u"));
        EXPECT_EQ(field_v(1,2,3), 4);
    }

    SECTION("test registry") {
        {
            Field field = MultiField {eckit::YAMLConfiguration{json()}}.field("temperature");
            auto temp = array::make_view<Value, 3>(field);
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
