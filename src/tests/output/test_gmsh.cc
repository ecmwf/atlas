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
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/output/Output.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"
#include "tests/TestMeshes.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

namespace {

struct NodeData {
    int components;
    std::vector<gidx_t> indices;
    std::vector<double> values;
};

NodeData read_scalar_node_data(const std::string& path) {
    std::ifstream file(path);
    std::string line;
    while (std::getline(file, line) && line != "$NodeData") {
    }
    EXPECT(file);

    int nb_string_tags;
    file >> nb_string_tags;
    std::getline(file, line);
    for (int tag = 0; tag < nb_string_tags; ++tag) {
        std::getline(file, line);
    }

    int nb_real_tags;
    file >> nb_real_tags;
    for (int tag = 0; tag < nb_real_tags; ++tag) {
        double value;
        file >> value;
    }

    int nb_integer_tags;
    file >> nb_integer_tags;
    std::vector<int> integer_tags(nb_integer_tags);
    for (int& tag : integer_tags) {
        file >> tag;
    }
    EXPECT_EQ(integer_tags.size(), 4);

    NodeData data;
    data.components      = integer_tags[1];
    const int nb_entries = integer_tags[2];
    data.indices.reserve(nb_entries);
    data.values.reserve(nb_entries * data.components);
    std::getline(file, line);
    for (int entry = 0; entry < nb_entries; ++entry) {
        std::getline(file, line);
        std::istringstream record(line);
        gidx_t index;
        record >> index;
        data.indices.push_back(index);
        for (int component = 0; component < data.components; ++component) {
            std::string value;
            record >> value;
            data.values.push_back(std::stod(value));
        }
    }
    return data;
}

Field field_with_missing_value() {
    Field field("values", array::make_datatype<double>(), array::make_shape(3));
    auto values = array::make_view<double, 1>(field);
    values(0)   = 1.;
    values(1)   = -999.;
    values(2)   = 3.;
    field.metadata().set("missing_value_type", "equals");
    field.metadata().set("missing_value", -999.);
    return field;
}

Field vector_field_with_missing_value() {
    Field field("vector_values", array::make_datatype<double>(), array::make_shape(2, 2));
    field.set_variables(2);
    auto values  = array::make_view<double, 2>(field);
    values(0, 0) = 1.;
    values(0, 1) = 2.;
    values(1, 0) = 3.;
    values(1, 1) = -999.;
    field.metadata().set("missing_value_type", "equals");
    field.metadata().set("missing_value", -999.);
    return field;
}

Field integer_field_with_missing_value() {
    Field field("integer_values", array::make_datatype<int>(), array::make_shape(2));
    auto values = array::make_view<int, 1>(field);
    values(0)   = 1;
    values(1)   = -999;
    field.metadata().set("missing_value_type", "equals");
    field.metadata().set("missing_value", -999);
    return field;
}

util::Config missing_value_config(const std::string& policy, double fill = 0.) {
    util::Config missing_value;
    missing_value.set("policy", policy);
    if (policy == "fill") {
        missing_value.set("fill", fill);
    }
    return util::Config("missing_value", missing_value);
}

}  // namespace

//-----------------------------------------------------------------------------

CASE("test_gmsh_output_1") {
    Mesh mesh = test::generate_mesh(Grid("N32"));
    output::Gmsh gmsh("test_gmsh_output_1.msh");
    gmsh.write(mesh);
}

CASE("test_gmsh_output_2") {
    Mesh mesh = test::generate_mesh(Grid("N32"));
    atlas::output::GmshFileStream file("test_gmsh_output_2.msh", "w");
    Log::warning() << "TODO: Not yet implemented!!! ATLAS-254" << std::endl;
    // output::Gmsh gmsh( file );
    // gmsh.write( mesh );
}

CASE("test_gmsh_missing_value_policy") {
    SECTION("skip is the default") {
        output::Gmsh("test_gmsh_missing_skip.msh").write(field_with_missing_value());
        const auto data = read_scalar_node_data("test_gmsh_missing_skip.msh");

        EXPECT_EQ(data.components, 1);
        EXPECT_EQ(data.indices, std::vector<gidx_t>({1, 3}));
        EXPECT_EQ(data.values, std::vector<double>({1., 3.}));
    }

    SECTION("fill") {
        output::Gmsh("test_gmsh_missing_fill.msh", missing_value_config("fill", -7.5))
            .write(field_with_missing_value());
        const auto data = read_scalar_node_data("test_gmsh_missing_fill.msh");

        EXPECT_EQ(data.components, 1);
        EXPECT_EQ(data.indices, std::vector<gidx_t>({1, 2, 3}));
        EXPECT_EQ(data.values, std::vector<double>({1., -7.5, 3.}));
    }

    SECTION("fill defaults to zero") {
        output::Gmsh("test_gmsh_missing_fill_default.msh", missing_value_config("fill"))
            .write(field_with_missing_value());
        const auto data = read_scalar_node_data("test_gmsh_missing_fill_default.msh");

        EXPECT_EQ(data.values, std::vector<double>({1., 0., 3.}));
    }

    SECTION("write config overrides constructor fill policy") {
        output::Gmsh gmsh("test_gmsh_missing_write_skip.msh", missing_value_config("fill", -7.5));
        gmsh.write(field_with_missing_value(), missing_value_config("skip"));
        const auto data = read_scalar_node_data("test_gmsh_missing_write_skip.msh");

        EXPECT_EQ(data.indices, std::vector<gidx_t>({1, 3}));
        EXPECT_EQ(data.values, std::vector<double>({1., 3.}));
    }

    SECTION("write config overrides constructor skip policy") {
        output::Gmsh gmsh("test_gmsh_missing_write_fill.msh", missing_value_config("skip"));
        gmsh.write(field_with_missing_value(), missing_value_config("fill", -7.5));
        const auto data = read_scalar_node_data("test_gmsh_missing_write_fill.msh");

        EXPECT_EQ(data.indices, std::vector<gidx_t>({1, 2, 3}));
        EXPECT_EQ(data.values, std::vector<double>({1., -7.5, 3.}));
    }

    SECTION("nan") {
        output::Gmsh("test_gmsh_missing_nan.msh", missing_value_config("nan")).write(field_with_missing_value());
        const auto data = read_scalar_node_data("test_gmsh_missing_nan.msh");

        EXPECT_EQ(data.components, 1);
        EXPECT_EQ(data.indices, std::vector<gidx_t>({1, 2, 3}));
        EXPECT_EQUAL(data.values[0], 1.);
        EXPECT(std::isnan(data.values[1]));
        EXPECT_EQUAL(data.values[2], 3.);
    }

    SECTION("preserve") {
        output::Gmsh("test_gmsh_missing_preserve.msh", missing_value_config("preserve"))
            .write(field_with_missing_value());
        const auto data = read_scalar_node_data("test_gmsh_missing_preserve.msh");

        EXPECT_EQ(data.indices, std::vector<gidx_t>({1, 2, 3}));
        EXPECT_EQ(data.values, std::vector<double>({1., -999., 3.}));
    }

    SECTION("fill replaces only missing vector components") {
        output::Gmsh("test_gmsh_missing_vector_fill.msh", missing_value_config("fill", -7.5))
            .write(vector_field_with_missing_value());
        const auto data = read_scalar_node_data("test_gmsh_missing_vector_fill.msh");

        EXPECT_EQ(data.components, 3);
        EXPECT_EQ(data.indices, std::vector<gidx_t>({1, 2}));
        EXPECT_EQ(data.values, std::vector<double>({1., 2., 0., 3., -7.5, 0.}));
    }

    SECTION("invalid policy is rejected") {
        EXPECT_THROWS(output::Gmsh("test_gmsh_missing_invalid.msh", missing_value_config("invalid"))
                          .write(field_with_missing_value()));
    }

    SECTION("nan is rejected for integer fields") {
        EXPECT_THROWS(output::Gmsh("test_gmsh_missing_integer_nan.msh", missing_value_config("nan"))
                          .write(integer_field_with_missing_value()));
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
