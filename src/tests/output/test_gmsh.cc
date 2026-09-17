/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/output/Gmsh.h"
#include "atlas/output/Output.h"
#include "atlas/util/Topology.h"
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

std::vector<gidx_t> read_element_tags(const std::string& path, int requested_dimension) {
    std::ifstream file(path);
    std::string line;
    while (std::getline(file, line) && line != "$Elements") {
    }
    EXPECT(file);

    size_t nb_blocks;
    size_t nb_elements;
    gidx_t min_tag;
    gidx_t max_tag;
    file >> nb_blocks >> nb_elements >> min_tag >> max_tag;
    std::getline(file, line);

    std::vector<gidx_t> tags;
    for (size_t block = 0; block < nb_blocks; ++block) {
        int dimension;
        int entity_tag;
        int element_type;
        size_t nb_elements_in_block;
        file >> dimension >> entity_tag >> element_type >> nb_elements_in_block;
        std::getline(file, line);
        for (size_t element = 0; element < nb_elements_in_block; ++element) {
            std::getline(file, line);
            if (dimension == requested_dimension) {
                std::istringstream record(line);
                gidx_t tag;
                record >> tag;
                tags.push_back(tag);
            }
        }
    }
    return tags;
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

util::Config masked_value_config(const std::string& policy, double fill = 0.) {
    util::Config masked_value;
    masked_value.set("policy", policy);
    if (policy == "fill") {
        masked_value.set("fill", fill);
    }
    return util::Config("masked_value", masked_value);
}

}  // namespace

//-----------------------------------------------------------------------------

CASE("test_gmsh_output_1") {
    Mesh mesh = test::generate_mesh(Grid("N32"));
    mesh.metadata().set("part", 2);
    mesh.metadata().set("nb_parts", 4);
    output::Gmsh gmsh("test_gmsh_output_1.msh", util::Config("element_partition_as_entity", true));
    gmsh.write(mesh);

    std::ifstream file("test_gmsh_output_1.msh");
    std::string line;
    std::getline(file, line);
    EXPECT_EQ(line, "$MeshFormat");
    std::getline(file, line);
    EXPECT_EQ(line, "4.1 0 " + std::to_string(sizeof(size_t)));

    while (std::getline(file, line) && line != "$Entities") {
    }
    EXPECT(file);
    int nb_points;
    int nb_curves;
    int nb_surfaces;
    int nb_volumes;
    file >> nb_points >> nb_curves >> nb_surfaces >> nb_volumes;
    EXPECT_EQ(nb_points, 0);
    EXPECT_EQ(nb_curves, 0);
    EXPECT_EQ(nb_surfaces, 2);
    EXPECT_EQ(nb_volumes, 0);

    std::getline(file, line);
    for (int expected_owner : {0, 2}) {
        std::getline(file, line);
        std::istringstream entity_record(line);
        int entity_tag;
        entity_record >> entity_tag;
        EXPECT_EQ(entity_tag, expected_owner);
    }

    while (std::getline(file, line) && line != "$Nodes") {
    }
    EXPECT(file);
    size_t nb_node_blocks;
    size_t nb_nodes;
    gidx_t min_node_tag;
    gidx_t max_node_tag;
    file >> nb_node_blocks >> nb_nodes >> min_node_tag >> max_node_tag;
    EXPECT_EQ(nb_node_blocks, 1);
    EXPECT_EQ(nb_nodes, mesh.nodes().size());

    int entity_dimension;
    int entity_tag;
    int parametric;
    size_t nb_nodes_in_block;
    file >> entity_dimension >> entity_tag >> parametric >> nb_nodes_in_block;
    EXPECT_EQ(entity_dimension, 2);
    EXPECT_EQ(entity_tag, 2);
    EXPECT_EQ(parametric, 0);
    EXPECT_EQ(nb_nodes_in_block, nb_nodes);
}

CASE("test_gmsh_output_binary") {
    Mesh mesh = test::generate_mesh(Grid("N16"));
    util::Config config("binary", true);
    config.set("info", true);
    output::Gmsh("test_gmsh_output_binary.msh", config).write(mesh);

    std::ifstream file("test_gmsh_output_binary.msh", std::ios::binary);
    std::string line;
    std::getline(file, line);
    EXPECT_EQ(line, "$MeshFormat");
    std::getline(file, line);
    EXPECT_EQ(line, "4.1 1 " + std::to_string(sizeof(size_t)));

    std::ifstream info_file("test_gmsh_output_binary_info.msh", std::ios::binary);
    std::getline(info_file, line);
    EXPECT_EQ(line, "$MeshFormat");
    std::getline(info_file, line);
    EXPECT_EQ(line, "4.1 1 " + std::to_string(sizeof(size_t)));
    int endian_check;
    info_file.read(reinterpret_cast<char*>(&endian_check), sizeof(endian_check));
    EXPECT_EQ(endian_check, 1);
    while (std::getline(info_file, line) && line != "$NodeData") {
    }
    EXPECT(info_file);
    for (int header_line = 0; header_line < 9; ++header_line) {
        std::getline(info_file, line);
    }
    int node_tag;
    double latitude;
    info_file.read(reinterpret_cast<char*>(&node_tag), sizeof(node_tag));
    info_file.read(reinterpret_cast<char*>(&latitude), sizeof(latitude));
    auto global_index = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
    auto lonlat       = array::make_view<double, 2>(mesh.nodes().lonlat());
    EXPECT_EQ(node_tag, global_index(0));
    EXPECT_EQ(latitude, lonlat(0, LAT));
}

CASE("test_gmsh_preserves_32bit_edge_tags") {
    Mesh mesh = test::generate_mesh(Grid("N16"));
    mesh::actions::build_edges(mesh);

    auto edge_global_index = array::make_view<gidx_t, 1>(mesh.edges().global_index());
    for (idx_t edge = 0; edge < mesh.edges().size(); ++edge) {
        edge_global_index(edge) = 10 + 1000 * edge;
    }

    output::Gmsh("test_gmsh_output_compact_edges.msh").write(mesh, util::Config("elements", "edges"));
    auto tags = read_element_tags("test_gmsh_output_compact_edges.msh", 1);
    std::sort(tags.begin(), tags.end());
    EXPECT(tags.size() >= 2);
    if (tags.size() < 2) {
        return;
    }

    EXPECT_EQ(tags[0], 1);
    EXPECT_EQ(tags[1], 1001);
}

CASE("test_gmsh_selects_cell_or_edge_elements") {
    Mesh mesh = test::generate_mesh(Grid("N16"));
    mesh::actions::build_edges(mesh);

    output::Gmsh("test_gmsh_output_auto.msh").write(mesh);
    EXPECT(!read_element_tags("test_gmsh_output_auto.msh", 2).empty());
    EXPECT(read_element_tags("test_gmsh_output_auto.msh", 1).empty());

    output::Gmsh("test_gmsh_output_cells.msh").write(mesh, util::Config("elements", "cells"));
    EXPECT(!read_element_tags("test_gmsh_output_cells.msh", 2).empty());
    EXPECT(read_element_tags("test_gmsh_output_cells.msh", 1).empty());

    output::Gmsh("test_gmsh_output_edges.msh").write(mesh, util::Config("elements", "edges"));
    EXPECT(read_element_tags("test_gmsh_output_edges.msh", 2).empty());
    EXPECT(!read_element_tags("test_gmsh_output_edges.msh", 1).empty());

    mesh.cells().clear();
    output::Gmsh("test_gmsh_output_auto_edges.msh").write(mesh);
    EXPECT(read_element_tags("test_gmsh_output_auto_edges.msh", 2).empty());
    EXPECT(!read_element_tags("test_gmsh_output_auto_edges.msh", 1).empty());

    mesh.edges().clear();
    output::Gmsh("test_gmsh_output_no_elements.msh").write(mesh);
    EXPECT(read_element_tags("test_gmsh_output_no_elements.msh", 2).empty());
    EXPECT(read_element_tags("test_gmsh_output_no_elements.msh", 1).empty());

    EXPECT_THROWS(output::Gmsh("test_gmsh_output_invalid.msh").write(mesh, util::Config("elements", "both")));
}

CASE("test_gmsh_output_canonical_entities") {
    Mesh mesh = test::generate_mesh(Grid("N16"));
    mesh.metadata().set("part", 2);
    mesh.metadata().set("nb_parts", 4);
    output::Gmsh("test_gmsh_output_canonical.msh", util::Config("element_partition_as_entity", false)).write(mesh);

    std::ifstream file("test_gmsh_output_canonical.msh");
    std::string contents((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    EXPECT(contents.find("$PartitionedEntities\n4\n0\n0 0 1 0\n4 2 1 1 2 ") != std::string::npos);
    EXPECT(contents.find("$Nodes\n1 ") != std::string::npos);
    EXPECT(contents.find("\n2 4 0 ") != std::string::npos);
}

CASE("test_gmsh_output_ghost_elements") {
    Mesh mesh = test::generate_mesh(Grid("N16"));
    mesh.metadata().set("part", 2);
    mesh.metadata().set("nb_parts", 4);

    mesh::Elements& elements = mesh.cells().elements(0);
    auto flags               = elements.view<int, 1>(elements.flags());
    auto halo                = elements.view<int, 1>(elements.halo());
    auto partition           = elements.view<int, 1>(elements.partition());
    util::Topology::set(flags(0), util::Topology::GHOST);
    halo(0)      = 1;
    partition(0) = 0;

    util::Config owner_entity_config("element_partition_as_entity", true);
    output::Gmsh("test_gmsh_output_without_ghost.msh", owner_entity_config).write(mesh);
    std::ifstream without_ghost_file("test_gmsh_output_without_ghost.msh");
    std::string without_ghost((std::istreambuf_iterator<char>(without_ghost_file)), std::istreambuf_iterator<char>());

    owner_entity_config.set("ghost", true);
    output::Gmsh("test_gmsh_output_ghost.msh", owner_entity_config).write(mesh);

    std::ifstream file("test_gmsh_output_ghost.msh");
    std::string contents((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    auto element_count = [](const std::string& mesh_contents) {
        std::istringstream stream(mesh_contents);
        std::string line;
        while (std::getline(stream, line) && line != "$Elements") {
        }
        size_t nb_blocks;
        size_t nb_elements;
        stream >> nb_blocks >> nb_elements;
        return nb_elements;
    };
    EXPECT_EQ(element_count(contents), element_count(without_ghost) + 1);
    EXPECT(without_ghost.find("$GhostElements") == std::string::npos);
    EXPECT(contents.find("$PartitionedEntities") == std::string::npos);

    const auto elements_begin = contents.find("$Elements\n");
    const auto elements_end   = contents.find("$EndElements\n", elements_begin);
    EXPECT(elements_begin != std::string::npos);
    EXPECT(elements_end != std::string::npos);
    EXPECT(contents.substr(elements_begin, elements_end - elements_begin).find("\n2 0 ") != std::string::npos);

    EXPECT(contents.find("$GhostElements") == std::string::npos);
    EXPECT(contents.find("$ElementData") == std::string::npos);

    util::Config canonical_config("ghost", true);
    canonical_config.set("element_partition_as_entity", false);
    output::Gmsh("test_gmsh_output_ghost_canonical.msh", canonical_config).write(mesh);

    std::ifstream canonical_file("test_gmsh_output_ghost_canonical.msh");
    std::string canonical((std::istreambuf_iterator<char>(canonical_file)), std::istreambuf_iterator<char>());
    EXPECT(canonical.find("$GhostElements") == std::string::npos);
    std::istringstream canonical_stream(canonical.substr(canonical.find("$Elements\n") + 10));
    size_t nb_blocks;
    size_t nb_elements;
    gidx_t min_element_tag;
    gidx_t max_element_tag;
    canonical_stream >> nb_blocks >> nb_elements >> min_element_tag >> max_element_tag;
    bool found_local_entity = false;
    bool found_owner_entity = false;
    for (size_t block = 0; block < nb_blocks; ++block) {
        int dimension;
        int entity_tag;
        int element_type;
        size_t block_size;
        canonical_stream >> dimension >> entity_tag >> element_type >> block_size;
        EXPECT_EQ(dimension, 2);
        found_local_entity |= entity_tag == 4;
        found_owner_entity |= entity_tag == 2;
        std::string line;
        std::getline(canonical_stream, line);
        for (size_t element = 0; element < block_size; ++element) {
            std::getline(canonical_stream, line);
        }
    }
    EXPECT(found_local_entity);
    EXPECT(found_owner_entity);

    util::Config binary_config("ghost", true);
    binary_config.set("binary", true);
    binary_config.set("element_partition_as_entity", true);
    output::Gmsh("test_gmsh_output_ghost_binary.msh", binary_config).write(mesh);

    std::ifstream binary_file("test_gmsh_output_ghost_binary.msh", std::ios::binary);
    std::getline(binary_file, contents);
    EXPECT_EQ(contents, "$MeshFormat");
    std::getline(binary_file, contents);
    EXPECT_EQ(contents, "4.1 1 " + std::to_string(sizeof(size_t)));
    std::string binary_contents((std::istreambuf_iterator<char>(binary_file)), std::istreambuf_iterator<char>());
    EXPECT(binary_contents.find("$ElementData") == std::string::npos);
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

CASE("test_gmsh_masked_value_policy") {
    auto write_masked_field = [](const std::string& path, const util::Config& config) {
        Mesh mesh = test::generate_mesh(Grid("O8"));
        functionspace::NodeColumns functionspace(mesh);
        functionspace.setMask([](int* mask, size_t size) {
            std::fill(mask, mask + size, 1);
            mask[1] = 0;
        });

        Field field = functionspace.createField<double>(option::name("values"));
        auto values = array::make_view<double, 1>(field);
        for (idx_t n = 0; n < values.size(); ++n) {
            values(n) = static_cast<double>(n + 1);
        }
        output::Gmsh(path).write(field, config);
        return read_scalar_node_data(path);
    };

    SECTION("skip is the default") {
        const auto data = write_masked_field("test_gmsh_masked_preserve.msh", util::Config());
        EXPECT(std::find(data.values.begin(), data.values.end(), 2.) == data.values.end());
    }

    SECTION("skip") {
        const auto data = write_masked_field("test_gmsh_masked_skip.msh", masked_value_config("skip"));
        EXPECT(std::find(data.values.begin(), data.values.end(), 2.) == data.values.end());
    }

    SECTION("fill") {
        const auto data = write_masked_field("test_gmsh_masked_fill.msh", masked_value_config("fill", -7.5));
        EXPECT_EQ(data.values[1], -7.5);
    }

    SECTION("fill takes precedence over missing-value skip") {
        Mesh mesh = test::generate_mesh(Grid("O8"));
        functionspace::NodeColumns functionspace(mesh);
        functionspace.setMask([](int* mask, size_t size) {
            std::fill(mask, mask + size, 1);
            mask[0] = 0;
        });
        Field field = functionspace.createField<double>(option::name("values"));
        auto values = array::make_view<double, 1>(field);
        std::fill(values.data(), values.data() + values.size(), 1.);
        values(0) = -999.;
        field.metadata().set("missing_value_type", "equals");
        field.metadata().set("missing_value", -999.);

        output::Gmsh("test_gmsh_masked_missing_fill.msh").write(field, masked_value_config("fill", -7.5));
        const auto data = read_scalar_node_data("test_gmsh_masked_missing_fill.msh");

        EXPECT_EQ(data.values[0], -7.5);
    }

    SECTION("masked and missing values use different fill values") {
        Mesh mesh = test::generate_mesh(Grid("O8"));
        functionspace::NodeColumns functionspace(mesh);
        functionspace.setMask([](int* mask, size_t size) {
            std::fill(mask, mask + size, 1);
            mask[0] = 0;
        });
        Field field = functionspace.createField<double>(option::name("values"));
        auto values = array::make_view<double, 1>(field);
        std::fill(values.data(), values.data() + values.size(), 1.);
        values(1)   = -999.;
        field.metadata().set("missing_value_type", "equals");
        field.metadata().set("missing_value", -999.);

        const util::Config config = masked_value_config("fill", -7.5) | missing_value_config("fill", -8.5);
        output::Gmsh("test_gmsh_masked_missing_different_fill.msh").write(field, config);
        const auto data = read_scalar_node_data("test_gmsh_masked_missing_different_fill.msh");

        EXPECT_EQ(data.values[0], -7.5);
        EXPECT_EQ(data.values[1], -8.5);
    }

    SECTION("nan") {
        const auto data = write_masked_field("test_gmsh_masked_nan.msh", masked_value_config("nan"));
        EXPECT(std::isnan(data.values[1]));
    }

    SECTION("invalid policy is rejected") {
        EXPECT_THROWS(write_masked_field("test_gmsh_masked_invalid.msh", masked_value_config("invalid")));
    }

    SECTION("nan is rejected for integer fields") {
        Mesh mesh = test::generate_mesh(Grid("O8"));
        functionspace::NodeColumns functionspace(mesh);
        functionspace.setMask([](int* mask, size_t size) {
            std::fill(mask, mask + size, 1);
            mask[0] = 0;
        });
        Field field = functionspace.createField<int>(option::name("values"));

        EXPECT_THROWS(output::Gmsh("test_gmsh_masked_integer_nan.msh")
                          .write(field, masked_value_config("nan")));
    }

    SECTION("does not allocate an unset mask") {
        Mesh mesh = test::generate_mesh(Grid("O8"));
        functionspace::NodeColumns functionspace(mesh);
        EXPECT(!functionspace.hasMask());

        Field field = functionspace.createField<double>(option::name("values"));
        auto values = array::make_view<double, 1>(field);
        std::fill(values.data(), values.data() + values.size(), 1.);
        output::Gmsh("test_gmsh_masked_unset.msh").write(field, masked_value_config("fill"));

        EXPECT(!functionspace.hasMask());
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
