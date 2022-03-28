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
    long nfld   = 0;
    long nlev   = 0;
    long nblk   = 0;

    if (!p.get("ngptot", ngptot)) {
        throw_Exception("Could not find 'ngptot' in Parametrisation");
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
    nfld            = fields.size();
    size_t nvar_tot = 0;
    for (const auto& field_config : fields) {
        size_t nvar = 1;
        field_config.get("nvar", nvar);
        nvar_tot += nvar;
    }

    auto multiarray_spec = array::make_shape(nblk, nvar_tot, nlev, nproma);
    auto& multiarray     = multifield.allocate(datatype, multiarray_spec);

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
        multifield.fields_.emplace_back(field);
        ATLAS_ASSERT(multifield.field_index_.find(field.name()) == multifield.field_index_.end(),
                     "Field with name already exists!");
        multifield.field_index_[field.name()] = i;
        multiarray_field_idx += field_vars;
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
        p.set("datatype", array::make_datatype<double>().str());
        p.set("fields", [] {
            std::vector<util::Config> fields(5);
            fields[0].set("name", "temperature");

            fields[1].set("name", "pressure");

            fields[2].set("name", "density");

            fields[3].set("name", "clv");
            fields[3].set("nvar", 5);

            fields[4].set("name", "wind_u");
            return fields;
        }());

        // We can also translate parameters to a json:
        std::stringstream json;
        eckit::JSON js(json);
        js << p;
        return json.str();
    };


    // And we can create back parameters from json:
    Log::info() << "json = " << json() << std::endl;
    std::stringstream json_stream;
    json_stream << json();
    util::Config config(json_stream);

    MultiField multifield("MultiFieldCreatorIFS", config);

    const auto nblk = multifield.array().shape(0);
    const auto nfld = multifield.size();
    EXPECT_EQ(nfld, 5);

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

    auto temp   = array::make_view<double, 3>(multifield.field("temperature"));
    auto pres   = array::make_view<double, 3>(multifield.field("pressure"));
    auto dens   = array::make_view<double, 3>(multifield.field("density"));
    auto clv    = array::make_view<double, 4>(multifield.field("clv"));  // note rank 4
    auto wind_u = array::make_view<double, 3>(multifield.field("wind_u"));

    // or
    {
        auto temp   = array::make_view<double, 3>(multifield[0]);
        auto pres   = array::make_view<double, 3>(multifield[1]);
        auto dens   = array::make_view<double, 3>(multifield[2]);
        auto clv    = array::make_view<double, 4>(multifield[3]);  // note rank 4
        auto wind_u = array::make_view<double, 3>(multifield[4]);
    }

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


    // Underlying array of MultiField
    {
        auto multiarray = array::make_view<double, 4>(multifield);
        EXPECT_EQ(multiarray.stride(0), block_stride);
        EXPECT_EQ(multiarray.stride(1), field_stride);
        EXPECT_EQ(multiarray.stride(2), level_stride);
        EXPECT_EQ(multiarray.stride(3), nproma_stride);

        EXPECT_EQ(multiarray(1, 0, 2, 3), 4.);
        EXPECT_EQ(multiarray(5, 1, 6, 7), 8.);
        EXPECT_EQ(multiarray(9, 2, 10, 11), 12.);
        EXPECT_EQ(multiarray(13, 5, 14, 15), 16.);
        EXPECT_EQ(multiarray(17, 8, 18, 3), 19.);


        EXPECT_EQ(multiarray.size(), nblk * (nfld - 1 + 5) * nlev * nproma);
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
