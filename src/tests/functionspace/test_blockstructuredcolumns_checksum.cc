/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <vector>
#include <cmath>

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Checksum.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas;
using namespace atlas::functionspace;
using namespace atlas::test;

namespace {

idx_t lane_value_count(const Field& field) {
    idx_t count = 1;
    for (idx_t dim = 1; dim < field.rank() - 1; ++dim) {
        count *= field.shape(dim);
    }
    return count;
}

template <typename Value>
void fill_block_values(Value* data, const Field& field, idx_t dim, idx_t offset, idx_t block_size, Value& next_value) {
    if (dim == field.rank() - 1) {
        for (idx_t jlane = 0; jlane < block_size; ++jlane) {
            data[offset + jlane * field.stride(dim)] = next_value;
            next_value += Value{1};
        }
        for (idx_t jlane = block_size; jlane < field.shape(dim); ++jlane) {
            data[offset + jlane * field.stride(dim)] = Value{-1};
        }
        return;
    }

    for (idx_t index = 0; index < field.shape(dim); ++index) {
        fill_block_values(data, field, dim + 1, offset + index * field.stride(dim), block_size, next_value);
    }
}

template <typename Value>
void pack_lane_values(std::vector<Value>& lane_values, const Field& field, const Value* data, idx_t dim, idx_t offset) {
    if (dim == field.rank() - 1) {
        lane_values.push_back(data[offset]);
        return;
    }

    for (idx_t index = 0; index < field.shape(dim); ++index) {
        pack_lane_values(lane_values, field, data, dim + 1, offset + index * field.stride(dim));
    }
}

template <typename Value>
std::string fill_field_and_expected_checksum(const BlockStructuredColumns& fs, Field& field) {
    auto* data = field.array().host_data<Value>();

    std::vector<util::checksum_t> lane_checksums;
    lane_checksums.reserve(fs.size());
    std::vector<Value> lane_values;
    lane_values.reserve(lane_value_count(field));

    Value next_value{1};
    for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
        auto blk = fs.block(jblk);
        const idx_t block_offset = jblk * field.stride(0);
        fill_block_values(data, field, 1, block_offset, blk.size(), next_value); // increments next_value

        for (idx_t jrof = 0; jrof < blk.size(); ++jrof) {
            lane_values.clear();
            if (field.rank() == 4) {
                // lev-outer, var-inner to match blocked checksum iteration order
                const idx_t nvar = field.shape(1);
                const idx_t nlev = field.shape(2);
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                        lane_values.push_back(data[block_offset + jvar * field.stride(1) + jlev * field.stride(2) + jrof * field.stride(3)]);
                    }
                }
            }
            else {
                pack_lane_values(lane_values, field, data, 1, block_offset + jrof * field.stride(field.rank() - 1));
            }
            lane_checksums.push_back(util::checksum(lane_values.data(), lane_values.size()));
        }
    }

    return util::checksum_to_hex_str(util::checksum(lane_checksums.data(), lane_checksums.size()));
}

template <typename Value>
std::string expected_fieldset_checksum(const BlockStructuredColumns& fs, const FieldSet& fields) {
    std::vector<util::checksum_t> lane_states(static_cast<size_t>(fs.size()));
    for (auto& lane_state : lane_states) {
        util::checksum_reset(lane_state);
    }

    std::vector<Value> lane_values;
    for (idx_t jfld = 0; jfld < fields.size(); ++jfld) {
        const Field& field = fields[jfld];
        const auto* data   = field.array().host_data<Value>();

        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            auto blk = fs.block(jblk);
            const idx_t block_offset = jblk * field.stride(0);

            for (idx_t jrof = 0; jrof < blk.size(); ++jrof) {
                lane_values.clear();
                if (field.rank() == 4) {
                    // lev-outer, var-inner to match blocked checksum iteration order
                    const idx_t nvar = field.shape(1);
                    const idx_t nlev = field.shape(2);
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                            lane_values.push_back(data[block_offset + jvar * field.stride(1) + jlev * field.stride(2) + jrof * field.stride(3)]);
                        }
                    }
                }
                else {
                    pack_lane_values(lane_values, field, data, 1, block_offset + jrof * field.stride(field.rank() - 1));
                }
                util::checksum_update(lane_states[static_cast<size_t>(blk.index(jrof))], lane_values.data(),
                                      lane_values.size());
            }
        }
    }

    std::vector<util::checksum_t> lane_digests(lane_states.size());
    for (size_t jlane = 0; jlane < lane_states.size(); ++jlane) {
        lane_digests[jlane] = util::checksum_digest(lane_states[jlane]);
    }
    return util::checksum_to_hex_str(util::checksum(lane_digests.data(), lane_digests.size()));
}

template <typename Value>
Field make_non_contiguous_field(const BlockStructuredColumns& fs, const util::Config& options, std::vector<Value>& rawdata) {
    Field prototype = fs.createField<Value>(option::name("field") | options);

    idx_t block_data_size = 1;
    for (idx_t dim = 1; dim < prototype.rank(); ++dim) {
        block_data_size *= prototype.shape(dim);
    }
    const idx_t block_stride = block_data_size + 5;

    array::ArrayShape shape;
    array::ArrayStrides strides;
    shape.reserve(prototype.rank());
    strides.resize(prototype.rank());
    for (idx_t dim = 0; dim < prototype.rank(); ++dim) {
        shape.emplace_back(prototype.shape(dim));
    }

    strides[prototype.rank() - 1] = 1;
    for (idx_t dim = prototype.rank() - 2; dim > 0; --dim) {
        strides[dim] = strides[dim + 1] * prototype.shape(dim + 1);
    }
    strides[0] = block_stride;

    rawdata.assign(fs.nblks() * block_stride, Value{-1});
    return Field("field", rawdata.data(), array::ArraySpec{shape, strides});
}

template <typename Value>
Field make_block_slice_non_contiguous_field(const BlockStructuredColumns& fs, const util::Config& options,
                                            std::vector<Value>& rawdata) {
    constexpr idx_t inner_dim_padding = 1;
    constexpr idx_t block_padding = 3;

    Field prototype = fs.createField<Value>(option::name("field") | options);

    array::ArrayShape shape = prototype.shape();
    array::ArrayStrides strides;
    strides.resize(prototype.rank());

    strides[prototype.rank() - 1] = 1;
    ATLAS_ASSERT(prototype.rank() > 2, "We always expect nproma to be contiguous");
    for (idx_t dim = prototype.rank() - 2; dim > 0; --dim) {
        // Add padding at each inner dimension so nproma-adjacent slices and higher dimensions are non-contiguous.
        strides[dim] = strides[dim + 1] * shape[dim + 1];
        strides[dim] += inner_dim_padding;
    }
    // Extra spacing between consecutive blocks in dimension 0.
    strides[0] = strides[1] * shape[1] + block_padding;

    rawdata.assign(fs.nblks() * strides[0], Value{-1});
    return Field("field", rawdata.data(), array::ArraySpec{shape, strides});
}

template <typename Value>
void expect_checksum_matches(const BlockStructuredColumns& fs, Field& field) {
    auto expected_checksum = fill_field_and_expected_checksum<Value>(fs, field);
    Log::info() << "Expected checksum: " << expected_checksum << std::endl;
    EXPECT_EQ(fs.checksum(field), expected_checksum);
}

template <typename Value>
void update_lane_states_with_global_field(const Field& field, std::vector<util::checksum_t>& lane_states) {
    ATLAS_ASSERT(field.shape(0) == static_cast<idx_t>(lane_states.size()));

    switch (field.rank()) {
        case 3: {
            auto value = array::make_view<Value, 3>(field);
            for (idx_t point = 0; point < field.shape(0); ++point) {
                util::checksum_update(lane_states[static_cast<size_t>(point)], &value(point, 0, 0), field.shape(1) * field.shape(2));
            }
            break;
        }
        case 2: {
            auto value = array::make_view<Value, 2>(field);
            for (idx_t point = 0; point < field.shape(0); ++point) {
                util::checksum_update(lane_states[static_cast<size_t>(point)], &value(point, 0), field.shape(1));
            }
            break;
        }
        case 1: {
            auto value = array::make_view<Value, 1>(field);
            for (idx_t point = 0; point < field.shape(0); ++point) {
                util::checksum_update(lane_states[static_cast<size_t>(point)], &value(point), 1);
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
util::checksum_t expected_global_field_checksum(const Field& global_field) {
    std::vector<util::checksum_t> lane_states(static_cast<size_t>(global_field.shape(0)));
    for (auto& lane_state : lane_states) {
        util::checksum_reset(lane_state);
    }

    update_lane_states_with_global_field<Value>(global_field, lane_states);

    return checksum_from_lane_states(lane_states);
}

template <typename Value>
util::checksum_t expected_global_fieldset_checksum(const FieldSet& global_fields) {
    ATLAS_ASSERT(global_fields.size() > 0);
    const idx_t npts = global_fields[0].shape(0);

    std::vector<util::checksum_t> lane_states(static_cast<size_t>(npts));
    for (auto& lane_state : lane_states) {
        util::checksum_reset(lane_state);
    }

    for (idx_t jfld = 0; jfld < global_fields.size(); ++jfld) {
        const Field& global = global_fields[jfld];
        update_lane_states_with_global_field<Value>(global, lane_states);
    }

    return checksum_from_lane_states(lane_states);
}

template <typename Value>
util::checksum_t fill_global_field_and_expected_checksum(Field& global) {
    auto* data = global.array().host_data<Value>();
    Value next_value{1};
    for (idx_t i = 0; i < global.size(); ++i) {
        data[i] = next_value;
        next_value += 1.;
    }
    return expected_global_field_checksum<Value>(global);
}

}  // namespace

namespace atlas {
namespace test {

CASE("test_BlockStructuredColumns checksum includes partial last block (mpi-serial)") {
    using Value = double;
    mpi::Scope scope("self"); // MPI-serial execution in this scope

    constexpr idx_t nlev = 4;
    constexpr idx_t nvar = 3;

    auto fs = [] {
        auto grid = StructuredGrid("O8");
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", grid.size() - 1);
        return BlockStructuredColumns(grid, config);
    }();
    REQUIRE(fs.size() > 1);
    REQUIRE(fs.size() % fs.nproma() == 1);
    REQUIRE(fs.block(fs.nblks() - 1).size() == 1);

    auto run_subcase = [&](const char* label, bool use_non_contiguous) {
        ATLAS_TRACE(label);
        Log::info() << "Running subcase: " << label << std::endl;
        if (use_non_contiguous) {
            std::vector<Value> rawdata;
            Field field = make_non_contiguous_field(fs, option::variables(nvar) | option::levels(nlev), rawdata);

            EXPECT(not field.contiguous());
            expect_checksum_matches<Value>(fs, field);
        }
        else {
            Field field = fs.createField<Value>(option::name("field") | option::variables(nvar) | option::levels(nlev));

            EXPECT(field.contiguous());
            expect_checksum_matches<Value>(fs, field);
        }
    };

    run_subcase("contiguous field", false);
    run_subcase("non-contiguous field", true);
}

CASE("test_BlockStructuredColumns checksum for various rank fields (mpi-serial)") {
    using Value = double;
    mpi::Scope scope("self");

    auto fs = [] {
        StructuredGrid grid{"O8"};
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", grid.size() - 1);
        return BlockStructuredColumns(grid, config);
    }();

    constexpr idx_t nvar = 3;
    constexpr idx_t nlev = 4;

    auto run_subcase = [&](const char* label, const util::Config& field_options, idx_t expected_rank) {
        ATLAS_TRACE(label);
        Log::info() << "Running variant: " << label << std::endl;

        Field contiguous = fs.createField<Value>(option::name("field") | field_options);
        EXPECT_EQ(contiguous.rank(), expected_rank);
        EXPECT(contiguous.contiguous());
        expect_checksum_matches<Value>(fs, contiguous);

        std::vector<Value> rawdata;
        Field non_contiguous = make_non_contiguous_field(fs, field_options, rawdata);
        EXPECT_EQ(non_contiguous.rank(), expected_rank);
        EXPECT(not non_contiguous.contiguous());
        expect_checksum_matches<Value>(fs, non_contiguous);
    };

    run_subcase("four-dimensional field with levels and variables", option::variables(nvar) | option::levels(nlev), 4);
    run_subcase("three-dimensional field without levels", option::variables(nvar), 3);
    run_subcase("three-dimensional field without variables", option::levels(nlev), 3);
    run_subcase("two-dimensional field without levels and variables", util::Config(), 2);
}

CASE("test_BlockStructuredColumns checksum in mpi-serial uses generic path for block-slice non-contiguous fields") {
    using Value = double;
    mpi::Scope scope("self");

    constexpr idx_t root = 0;
    constexpr idx_t nlev = 4;
    constexpr idx_t nvar = 3;

    auto fs = [] {
        auto grid = StructuredGrid("O8");
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", 5);
        return BlockStructuredColumns(grid, config);
    }();

    Field global = fs.createField<Value>(option::name("global") | option::global(root) | option::variables(nvar) |
                                         option::levels(nlev));
    util::checksum_t expected_checksum = 0;
    if (mpi::comm().rank() == root) {
        expected_checksum = fill_global_field_and_expected_checksum<Value>(global);
    }
    mpi::comm().broadcast(expected_checksum, root);

    std::vector<Value> rawdata;
    Field local = make_block_slice_non_contiguous_field(fs, option::variables(nvar) | option::levels(nlev), rawdata);
    EXPECT(not local.contiguous());
    EXPECT(local.stride(local.rank() - 2) != local.shape(local.rank() - 1));

    fs.scatter(global, local);

    auto expected_checksum_hex_str = util::checksum_to_hex_str(expected_checksum);
    Log::info() << "Expected checksum: " << expected_checksum_hex_str << std::endl;
    EXPECT_EQ(fs.checksum(local), expected_checksum_hex_str);
}

CASE("test_BlockStructuredColumns checksum(fieldset) matches checksum(field) for a single field (mpi-serial)") {
    using Value = double;
    mpi::Scope scope("self");

    constexpr idx_t nlev = 4;
    constexpr idx_t nvar = 3;

    auto fs = [] {
        auto grid = StructuredGrid("O8");
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", grid.size() - 1);
        return BlockStructuredColumns(grid, config);
    }();

    {
        Field field = fs.createField<Value>(option::name("field") | option::variables(nvar) | option::levels(nlev));
        fill_field_and_expected_checksum<Value>(fs, field);

        FieldSet fieldset;
        fieldset.add(field);

        EXPECT_EQ(fs.checksum(fieldset), fs.checksum(field));
    }

    {
        std::vector<Value> rawdata;
        Field field = make_non_contiguous_field<Value>(fs, option::variables(nvar) | option::levels(nlev), rawdata);
        fill_field_and_expected_checksum<Value>(fs, field);

        FieldSet fieldset;
        fieldset.add(field);

        EXPECT_EQ(fs.checksum(fieldset), fs.checksum(field));
    }
}

CASE("test_BlockStructuredColumns checksum(fieldset) combines all fields as lane checksum updates (mpi-serial)") {
    using Value = double;
    mpi::Scope scope("self");

    constexpr idx_t nlev = 4;
    constexpr idx_t nvar = 3;

    auto fs = [] {
        auto grid = StructuredGrid("O8");
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", 5);
        return BlockStructuredColumns(grid, config);
    }();

    Field contiguous = fs.createField<Value>(option::name("contiguous") | option::variables(nvar) | option::levels(nlev));
    std::vector<Value> rawdata;
    Field non_contiguous = make_non_contiguous_field<Value>(fs, option::variables(nvar) | option::levels(nlev), rawdata);

    fill_field_and_expected_checksum<Value>(fs, contiguous);
    fill_field_and_expected_checksum<Value>(fs, non_contiguous);

    FieldSet fieldset;
    fieldset.add(contiguous);
    fieldset.add(non_contiguous);

    const std::string expected = expected_fieldset_checksum<Value>(fs, fieldset);
    Log::info() << "Expected fieldset checksum: " << expected << std::endl;
    EXPECT_EQ(fs.checksum(fieldset), expected);
}

template <typename Value>
void run_mpi_checksum_caseT(const functionspace::BlockStructuredColumns& fs, idx_t nlev, idx_t nvar) {
    constexpr idx_t root = 0;

    util::Config option_nlev_nvars;
    if (nvar > 0) {
        option_nlev_nvars.set(option::variables(nvar));
    }
    if (nlev > 0) {
        option_nlev_nvars.set(option::levels(nlev));
    }

    Field global = fs.createField<Value>(option::name("global") | option::global(root) | option_nlev_nvars);
    util::checksum_t expected_checksum = 0;
    if (mpi::comm().rank() == root) {
        expected_checksum = fill_global_field_and_expected_checksum<Value>(global);
    }
    mpi::comm().broadcast(expected_checksum, root);

    Field local = fs.createField<Value>(option::name("local") | option_nlev_nvars);
    fs.scatter(global, local);
    idx_t expected_rank = idx_t{2} /*nblks, nproma*/ + std::min<idx_t>(nlev,1) + std::min<idx_t>(nvar,1);
    EXPECT_EQ(local.rank(), expected_rank);

    auto expected_checksum_hex_str = util::checksum_to_hex_str(expected_checksum);
    Log::info() << "Expected checksum: " << expected_checksum_hex_str << std::endl;
    EXPECT_EQ(fs.checksum(local), expected_checksum_hex_str);
};

void run_mpi_checksum_cases(const functionspace::BlockStructuredColumns& fs, idx_t nlev, idx_t nvar) {
    run_mpi_checksum_caseT<double>(fs, nlev, nvar);
    run_mpi_checksum_caseT<float>(fs, nlev, nvar);
    run_mpi_checksum_caseT<int>(fs, nlev, nvar);
    run_mpi_checksum_caseT<long>(fs, nlev, nvar);
}

template <typename Value>
void run_mpi_fieldset_checksum_caseT(const functionspace::BlockStructuredColumns& fs, idx_t nlev, idx_t nvar) {
    constexpr idx_t root = 0;

    util::Config option_nlev_nvars;
    if (nvar > 0) {
        option_nlev_nvars.set(option::variables(nvar));
    }
    if (nlev > 0) {
        option_nlev_nvars.set(option::levels(nlev));
    }

    Field global_a = fs.createField<Value>(option::name("global_a") | option::global(root) | option_nlev_nvars);
    Field global_b = fs.createField<Value>(option::name("global_b") | option::global(root) | option_nlev_nvars);

    util::checksum_t expected_fieldset = 0;
    if (mpi::comm().rank() == root) {
        fill_global_field_and_expected_checksum<Value>(global_a);
        fill_global_field_and_expected_checksum<Value>(global_b);

        FieldSet global_fieldset;
        global_fieldset.add(global_a);
        global_fieldset.add(global_b);
        expected_fieldset = expected_global_fieldset_checksum<Value>(global_fieldset);
    }
    mpi::comm().broadcast(expected_fieldset, root);

    Field local_a = fs.createField<Value>(option::name("local_a") | option_nlev_nvars);
    Field local_b = fs.createField<Value>(option::name("local_b") | option_nlev_nvars);
    fs.scatter(global_a, local_a);
    fs.scatter(global_b, local_b);

    FieldSet local_fieldset;
    local_fieldset.add(local_a);
    local_fieldset.add(local_b);

    auto expected_fieldset_hex_str = util::checksum_to_hex_str(expected_fieldset);
    Log::info() << "Expected fieldset checksum (MPI): " << expected_fieldset_hex_str << std::endl;
    EXPECT_EQ(fs.checksum(local_fieldset), expected_fieldset_hex_str);
}

void run_mpi_fieldset_checksum_cases(const functionspace::BlockStructuredColumns& fs, idx_t nlev, idx_t nvar) {
    run_mpi_fieldset_checksum_caseT<double>(fs, nlev, nvar);
    run_mpi_fieldset_checksum_caseT<float>(fs, nlev, nvar);
    run_mpi_fieldset_checksum_caseT<int>(fs, nlev, nvar);
    run_mpi_fieldset_checksum_caseT<long>(fs, nlev, nvar);
}


CASE("test_BlockStructuredColumns checksum is deterministic in MPI") {
    auto fs = [] {
        auto grid = StructuredGrid("O8");
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", 8);
        return BlockStructuredColumns(grid, config);
    }();

    constexpr idx_t nlev = 4;
    constexpr idx_t nvar = 3;

    auto run_subcase = [&](const char* label, idx_t levels, idx_t vars) {
        Log::info() << "Running subcase: " << label << std::endl;
        run_mpi_checksum_cases(fs, levels, vars);
    };

    run_subcase("with variables and levels", nlev, nvar);
    run_subcase("with variables and no levels", 0, nvar);
    run_subcase("with levels and no variables", nlev, 0);
    run_subcase("with no variables and no levels", 0, 0);
}

CASE("test_BlockStructuredColumns checksum(fieldset) is deterministic in MPI") {
    auto fs = [] {
        auto grid = StructuredGrid("O8");
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", 8);
        return BlockStructuredColumns(grid, config);
    }();

    constexpr idx_t nlev = 4;
    constexpr idx_t nvar = 3;

    auto run_subcase = [&](const char* label, idx_t levels, idx_t vars) {
        Log::info() << "Running fieldset subcase: " << label << std::endl;
        run_mpi_fieldset_checksum_cases(fs, levels, vars);
    };

    run_subcase("with variables and levels", nlev, nvar);
    run_subcase("with variables and no levels", 0, nvar);
    run_subcase("with levels and no variables", nlev, 0);
    run_subcase("with no variables and no levels", 0, 0);
}

CASE("test_BlockStructuredColumns checksum with halo=2 (mpi-serial)") {
    using Value = double;
    mpi::Scope scope("self");

    constexpr idx_t nlev = 4;
    constexpr idx_t nvar = 3;

    auto grid = StructuredGrid("O8");
    auto fs_h0 = [&] {
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", 8);
        return BlockStructuredColumns(grid, config);
    }();
    auto fs_h2 = [&] {
        util::Config config;
        config.set("halo", 2);
        config.set("nproma", 8);
        return BlockStructuredColumns(grid, config);
    }();

    Field field_h0 = fs_h0.createField<Value>(option::name("field_h0") | option::variables(nvar) | option::levels(nlev));
    fill_field_and_expected_checksum<Value>(fs_h0, field_h0);

    Field global = fs_h0.createField<Value>(option::name("global") | option::global(0) | option::variables(nvar) | option::levels(nlev));
    fs_h0.gather(field_h0, global);

    Field field_h2 = fs_h2.createField<Value>(option::name("field_h2") | option::variables(nvar) | option::levels(nlev));
    fs_h2.scatter(global, field_h2);

    EXPECT_EQ(fs_h0.checksum(field_h0), fs_h2.checksum(field_h2));

    std::vector<Value> rawdata;
    Field non_contiguous = make_non_contiguous_field<Value>(fs_h2, option::variables(nvar) | option::levels(nlev), rawdata);
    EXPECT(not non_contiguous.contiguous());
    fill_field_and_expected_checksum<Value>(fs_h2, non_contiguous);
    EXPECT_EQ(fs_h2.checksum(non_contiguous), fs_h2.checksum(non_contiguous));
}

CASE("test_BlockStructuredColumns checksum(fieldset) with halo=2 (mpi-serial)") {
    using Value = double;
    mpi::Scope scope("self");

    constexpr idx_t nlev = 4;
    constexpr idx_t nvar = 3;

    auto grid = StructuredGrid("O8");
    auto fs_h0 = [&] {
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", 8);
        return BlockStructuredColumns(grid, config);
    }();
    auto fs_h2 = [&] {
        util::Config config;
        config.set("halo", 2);
        config.set("nproma", 8);
        return BlockStructuredColumns(grid, config);
    }();

    Field field_h0_a = fs_h0.createField<Value>(option::name("field_h0_a") | option::variables(nvar) | option::levels(nlev));
    Field field_h0_b = fs_h0.createField<Value>(option::name("field_h0_b") | option::variables(2));
    fill_field_and_expected_checksum<Value>(fs_h0, field_h0_a);
    fill_field_and_expected_checksum<Value>(fs_h0, field_h0_b);

    Field global_a = fs_h0.createField<Value>(option::name("global_a") | option::global(0) | option::variables(nvar) | option::levels(nlev));
    Field global_b = fs_h0.createField<Value>(option::name("global_b") | option::global(0) | option::variables(2));
    fs_h0.gather(field_h0_a, global_a);
    fs_h0.gather(field_h0_b, global_b);

    Field field_h2_a = fs_h2.createField<Value>(option::name("field_h2_a") | option::variables(nvar) | option::levels(nlev));
    Field field_h2_b = fs_h2.createField<Value>(option::name("field_h2_b") | option::variables(2));
    fs_h2.scatter(global_a, field_h2_a);
    fs_h2.scatter(global_b, field_h2_b);

    FieldSet fieldset_h0;
    FieldSet fieldset_h2;
    fieldset_h0.add(field_h0_a);
    fieldset_h0.add(field_h0_b);
    fieldset_h2.add(field_h2_a);
    fieldset_h2.add(field_h2_b);

    EXPECT_EQ(fs_h0.checksum(fieldset_h0), fs_h2.checksum(fieldset_h2));
}

CASE("test_BlockStructuredColumns checksum with halo=2 matches halo=0 exactly (mpi-serial)") {
    using Value = double;
    mpi::Scope scope("self");

    auto grid = StructuredGrid("O8");
    auto fs_h0 = [&] {
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", 8);
        return BlockStructuredColumns(grid, config);
    }();
    auto fs_h2 = [&] {
        util::Config config;
        config.set("halo", 2);
        config.set("nproma", 8);
        return BlockStructuredColumns(grid, config);
    }();

    Field field_h0 = fs_h0.createField<Value>(option::name("field_h0") | option::variables(3) | option::levels(4));
    fill_field_and_expected_checksum<Value>(fs_h0, field_h0);

    Field global = fs_h0.createField<Value>(option::name("global") | option::global(0) | option::variables(3) | option::levels(4));
    fs_h0.gather(field_h0, global);

    Field field_h2 = fs_h2.createField<Value>(option::name("field_h2") | option::variables(3) | option::levels(4));
    fs_h2.scatter(global, field_h2);
    EXPECT_EQ(fs_h0.checksum(field_h0), fs_h2.checksum(field_h2));

    Field field_h0_b = fs_h0.createField<Value>(option::name("field_h0_b") | option::variables(2));
    fill_field_and_expected_checksum<Value>(fs_h0, field_h0_b);
    Field global_b = fs_h0.createField<Value>(option::name("global_b") | option::global(0) | option::variables(2));
    fs_h0.gather(field_h0_b, global_b);
    Field field_h2_b = fs_h2.createField<Value>(option::name("field_h2_b") | option::variables(2));
    fs_h2.scatter(global_b, field_h2_b);

    FieldSet fieldset_h0;
    FieldSet fieldset_h2;
    fieldset_h0.add(field_h0);
    fieldset_h0.add(field_h0_b);
    fieldset_h2.add(field_h2);
    fieldset_h2.add(field_h2_b);
    EXPECT_EQ(fs_h0.checksum(fieldset_h0), fs_h2.checksum(fieldset_h2));
}

CASE("test_BlockStructuredColumns checksum with halo=2 matches halo=0 exactly (MPI)") {
    using Value = double;

    constexpr idx_t root = 0;
    auto grid = StructuredGrid("O8");
    auto fs_h0 = [&] {
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", 8);
        return BlockStructuredColumns(grid, config);
    }();
    auto fs_h2 = [&] {
        util::Config config;
        config.set("halo", 2);
        config.set("nproma", 8);
        return BlockStructuredColumns(grid, config);
    }();

    auto check_field_checksum_stable_with_halo = [&](const util::Config& options, Value first_value) {
        Field global = fs_h0.createField<Value>(option::name("global") | option::global(root) | options);
        if (mpi::comm().rank() == root) {
            fill_global_field_and_expected_checksum<Value>(global);
            // restart with requested first value for a deterministic variant per case
            auto* data = global.array().host_data<Value>();
            Value next = first_value;
            for (idx_t i = 0; i < global.size(); ++i) {
                data[i] = next;
                next += Value{1};
            }
        }

        Field local_h0 = fs_h0.createField<Value>(option::name("local_h0") | options);
        Field local_h2 = fs_h2.createField<Value>(option::name("local_h2") | options);
        fs_h0.scatter(global, local_h0);
        fs_h2.scatter(global, local_h2);
        EXPECT_EQ(fs_h0.checksum(local_h0), fs_h2.checksum(local_h2));
    };

    check_field_checksum_stable_with_halo(option::variables(3) | option::levels(4), Value{24000});
    check_field_checksum_stable_with_halo(option::variables(3), Value{25000});
    check_field_checksum_stable_with_halo(option::levels(4), Value{26000});
    check_field_checksum_stable_with_halo(util::Config{}, Value{27000});

    auto check_fieldset_checksum_stable_with_halo = [&](const util::Config& options_a, const util::Config& options_b) {
        Field global_a = fs_h0.createField<Value>(option::name("global_a") | option::global(root) |
                                                option::variables(2) | option::levels(3));
        Field global_b = fs_h0.createField<Value>(option::name("global_b") | option::global(root) | option::variables(4));
        if (mpi::comm().rank() == root) {
            fill_global_field_and_expected_checksum<Value>(global_a);
            fill_global_field_and_expected_checksum<Value>(global_b);
        }

        Field local_h0_a = fs_h0.createField<Value>(option::name("local_h0_a") | option::variables(2) | option::levels(3));
        Field local_h0_b = fs_h0.createField<Value>(option::name("local_h0_b") | option::variables(4));
        Field local_h2_a = fs_h2.createField<Value>(option::name("local_h2_a") | option::variables(2) | option::levels(3));
        Field local_h2_b = fs_h2.createField<Value>(option::name("local_h2_b") | option::variables(4));

        fs_h0.scatter(global_a, local_h0_a);
        fs_h0.scatter(global_b, local_h0_b);
        fs_h2.scatter(global_a, local_h2_a);
        fs_h2.scatter(global_b, local_h2_b);

        FieldSet fieldset_h0;
        FieldSet fieldset_h2;
        fieldset_h0.add(local_h0_a);
        fieldset_h0.add(local_h0_b);
        fieldset_h2.add(local_h2_a);
        fieldset_h2.add(local_h2_b);
        EXPECT_EQ(fs_h0.checksum(fieldset_h0), fs_h2.checksum(fieldset_h2));
    };

    check_fieldset_checksum_stable_with_halo(/*a*/ option::variables(2) | option::levels(3), /*b*/ option::variables(4));

}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}