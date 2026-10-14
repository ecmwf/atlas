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
#include <vector>

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Checksum.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas;
using namespace atlas::functionspace;
using namespace atlas::test;

namespace {

template <typename Value>
Value point_value(gidx_t global_index, idx_t dim1 = 0, idx_t dim2 = 0) {
    return static_cast<Value>(global_index * 100 + dim1 * 10 + dim2 + 1);
}

template <typename FunctionSpace, typename Value>
void fill_owned_and_ghost_values(const FunctionSpace& fs, Field& field, Value ghost_value) {
    auto global_index = array::make_view<gidx_t, 1>(fs.global_index());
    auto ghost        = array::make_view<int, 1>(fs.ghost());

    switch (field.rank()) {
        case 1: {
            auto view = array::make_view<Value, 1>(field);
            for (idx_t point = 0; point < field.shape(0); ++point) {
                view(point) = ghost(point) ? ghost_value : point_value<Value>(global_index(point));
            }
            break;
        }
        case 2: {
            auto view = array::make_view<Value, 2>(field);
            for (idx_t point = 0; point < field.shape(0); ++point) {
                for (idx_t dim1 = 0; dim1 < field.shape(1); ++dim1) {
                    view(point, dim1) = ghost(point) ? ghost_value : point_value<Value>(global_index(point), dim1);
                }
            }
            break;
        }
        case 3: {
            auto view = array::make_view<Value, 3>(field);
            for (idx_t point = 0; point < field.shape(0); ++point) {
                for (idx_t dim1 = 0; dim1 < field.shape(1); ++dim1) {
                    for (idx_t dim2 = 0; dim2 < field.shape(2); ++dim2) {
                        view(point, dim1, dim2) = ghost(point) ? ghost_value : point_value<Value>(global_index(point), dim1, dim2);
                    }
                }
            }
            break;
        }
        default:
            ATLAS_NOTIMPLEMENTED;
    }
    field.set_dirty();
}

template <typename Value>
void update_lane_states_with_global_field(const Field& field, std::vector<util::checksum_t>& lane_states) {
    ATLAS_ASSERT(field.shape(0) == static_cast<idx_t>(lane_states.size()));

    switch (field.rank()) {
        case 1: {
            auto view = array::make_view<Value, 1>(field);
            for (idx_t point = 0; point < field.shape(0); ++point) {
                util::checksum_update(lane_states[static_cast<size_t>(point)], &view(point), 1);
            }
            break;
        }
        case 2: {
            auto view = array::make_view<Value, 2>(field);
            for (idx_t point = 0; point < field.shape(0); ++point) {
                util::checksum_update(lane_states[static_cast<size_t>(point)], &view(point, 0),
                                      static_cast<size_t>(field.shape(1)));
            }
            break;
        }
        case 3: {
            auto view = array::make_view<Value, 3>(field);
            for (idx_t point = 0; point < field.shape(0); ++point) {
                util::checksum_update(lane_states[static_cast<size_t>(point)], &view(point, 0, 0),
                                      static_cast<size_t>(field.shape(1) * field.shape(2)));
            }
            break;
        }
        default:
            ATLAS_NOTIMPLEMENTED;
    }
}

util::checksum_t checksum_from_lane_states(const std::vector<util::checksum_t>& lane_states) {
    std::vector<util::checksum_t> lane_digests(lane_states.size());
    for (size_t point = 0; point < lane_states.size(); ++point) {
        lane_digests[point] = util::checksum_digest(lane_states[point]);
    }
    return util::checksum(lane_digests.data(), lane_digests.size());
}

template <typename Value>
std::string expected_field_checksum(const Field& global_field, idx_t root = 0) {
    util::checksum_t expected_checksum = 0;
    auto& comm = mpi::comm();

    if (comm.rank() == root) {
        std::vector<util::checksum_t> lane_states(static_cast<size_t>(global_field.shape(0)));
        for (auto& lane_state : lane_states) {
            util::checksum_reset(lane_state);
        }

        update_lane_states_with_global_field<Value>(global_field, lane_states);
        expected_checksum = checksum_from_lane_states(lane_states);
    }

    comm.broadcast(expected_checksum, root);
    return util::checksum_to_hex_str(expected_checksum);
}

template <typename Value>
std::string expected_fieldset_checksum(const FieldSet& global_fields, idx_t root = 0) {
    util::checksum_t expected_checksum = 0;
    auto& comm = mpi::comm();

    if (comm.rank() == root) {
        ATLAS_ASSERT(global_fields.size() > 0);
        std::vector<util::checksum_t> lane_states(static_cast<size_t>(global_fields[0].shape(0)));
        for (auto& lane_state : lane_states) {
            util::checksum_reset(lane_state);
        }

        for (idx_t field_idx = 0; field_idx < global_fields.size(); ++field_idx) {
            update_lane_states_with_global_field<Value>(global_fields[field_idx], lane_states);
        }

        expected_checksum = checksum_from_lane_states(lane_states);
    }

    comm.broadcast(expected_checksum, root);
    return util::checksum_to_hex_str(expected_checksum);
}

template <typename FunctionSpace>
Field gather_global_field(const FunctionSpace& fs, const Field& local_field, idx_t root = 0) {
    Field global_field = fs.createField(local_field, option::global(root));
    fs.gather(local_field, global_field);
    return global_field;
}

template <typename FunctionSpace>
FieldSet gather_global_fieldset(const FunctionSpace& fs, const FieldSet& local_fields, idx_t root = 0) {
    FieldSet global_fields;
    for (idx_t field_idx = 0; field_idx < local_fields.size(); ++field_idx) {
        global_fields.add(fs.createField(local_fields[field_idx], option::global(root)));
    }
    fs.gather(local_fields, global_fields);
    return global_fields;
}

Mesh generate_mesh() {
    Grid grid{"O16"};
    return StructuredMeshGenerator().generate(grid);
}

template <typename FunctionSpace>
void expect_checksum_matches_global_field(const char* expected_checksum) {
    auto levels = 3;
    FunctionSpace fs(generate_mesh(), option::halo(1) | option::levels(levels));
    Field field = fs.template createField<double>(option::name("field") | option::variables(2));

    auto first_ghost_value = -999.;
    fill_owned_and_ghost_values(fs, field, first_ghost_value);

    Field global_field = gather_global_field(fs, field);
    const auto checksum = expected_field_checksum<double>(global_field);
    EXPECT_EQ(fs.checksum(field), checksum);

    auto second_ghost_value = -12345.;
    fill_owned_and_ghost_values(fs, field, second_ghost_value);
    EXPECT_EQ(fs.checksum(field), checksum);
    EXPECT_EQ(checksum, expected_checksum);
}

template <typename FunctionSpace>
void expect_checksum_matches_global_fieldset(const char* expected_checksum) {
    auto levels = 2;
    FunctionSpace fs(generate_mesh(), option::halo(1) | option::levels(levels));
    Field field1 = fs.template createField<double>(option::name("field1"));
    Field field2 = fs.template createField<double>(option::name("field2") | option::variables(2));

    fill_owned_and_ghost_values(fs, field1, -1.);
    fill_owned_and_ghost_values(fs, field2, -2.);

    FieldSet fields;
    fields.add(field1);
    fields.add(field2);

    FieldSet global_fields = gather_global_fieldset(fs, fields);
    const auto checksum = expected_fieldset_checksum<double>(global_fields);
    EXPECT_EQ(fs.checksum(fields), checksum);
    EXPECT_EQ(checksum, expected_checksum);
}

}  // namespace

namespace atlas {
namespace test {

CASE("test_nodecolumns_checksum_matches_global_field") {
    expect_checksum_matches_global_field<NodeColumns>("7290");
}

CASE("test_nodecolumns_checksum_fieldset_matches_global_fieldset") {
    expect_checksum_matches_global_fieldset<NodeColumns>("b0c2");
}

CASE("test_cellcolumns_checksum_matches_global_field") {
    expect_checksum_matches_global_field<CellColumns>("f645");
}

CASE("test_cellcolumns_checksum_fieldset_matches_global_fieldset") {
    expect_checksum_matches_global_fieldset<CellColumns>("893a");
}

CASE("test_edgecolumns_checksum_matches_global_field") {
    expect_checksum_matches_global_field<EdgeColumns>("7098");
}

CASE("test_edgecolumns_checksum_fieldset_matches_global_fieldset") {
    expect_checksum_matches_global_fieldset<EdgeColumns>("263b");
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}