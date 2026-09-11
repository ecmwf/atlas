/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/util/relayout.h"

#include "eckit/config/Resource.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

namespace {

bool on_device() {
    static bool value = eckit::Resource<bool>("--on-device", false);
    return value;
}

void prepare_source(Field& field) {
    if (on_device()) {
        field.syncDevice();
    }
}

void prepare_target(Field& field) {
    if (on_device()) {
        field.allocateDevice();
    }
}

void prepare_target(FieldSet& fields) {
    if (on_device()) {
        fields.allocateDevice();
    }
}

void sync_target(Field& field) {
    if (on_device()) {
        field.syncHost();
    }
}

void sync_target(FieldSet& fields) {
    if (on_device()) {
        fields.syncHost();
    }
}

template <typename Value>
void fill_structured_field(Field& field, Value first_value) {
    Value next_value = first_value;
    if (field.variables() && field.levels()) {
        auto value = array::make_view<Value, 3>(field);
        for (idx_t point = 0; point < field.shape(0); ++point) {
            for (idx_t jlev = 0; jlev < field.shape(1); ++jlev) {
                for (idx_t jvar = 0; jvar < field.shape(2); ++jvar) {
                    value(point, jlev, jvar) = next_value;
                    next_value += Value{1};
                }
            }
        }
    }
    else if (field.variables() || field.levels()) {
        auto value = array::make_view<Value, 2>(field);
        for (idx_t point = 0; point < field.shape(0); ++point) {
            for (idx_t entry = 0; entry < field.shape(1); ++entry) {
                value(point, entry) = next_value;
                next_value += Value{1};
            }
        }
    }
    else {
        auto value = array::make_view<Value, 1>(field);
        for (idx_t point = 0; point < field.shape(0); ++point) {
            value(point) = next_value;
            next_value += Value{1};
        }
    }
}

template <typename Value>
void fill_blocked_field(const functionspace::BlockStructuredColumns& functionspace, Field& field, Value first_value) {
    Value next_value = first_value;
    if (field.variables() && field.levels()) {
        auto value = array::make_view<Value, 4>(field);
        for (idx_t jblk = 0; jblk < field.shape(0); ++jblk) {
            for (idx_t jlane = 0; jlane < functionspace.block(jblk).size(); ++jlane) {
                for (idx_t jlev = 0; jlev < field.shape(2); ++jlev) {
                    for (idx_t jvar = 0; jvar < field.shape(1); ++jvar) {
                        value(jblk, jvar, jlev, jlane) = next_value;
                        next_value += Value{1};
                    }
                }
            }
        }
    }
    else if (field.variables() || field.levels()) {
        auto value = array::make_view<Value, 3>(field);
        for (idx_t jblk = 0; jblk < field.shape(0); ++jblk) {
            for (idx_t jlane = 0; jlane < functionspace.block(jblk).size(); ++jlane) {
                for (idx_t entry = 0; entry < field.shape(1); ++entry) {
                    value(jblk, entry, jlane) = next_value;
                    next_value += Value{1};
                }
            }
        }
    }
    else {
        auto value = array::make_view<Value, 2>(field);
        for (idx_t jblk = 0; jblk < field.shape(0); ++jblk) {
            for (idx_t jlane = 0; jlane < functionspace.block(jblk).size(); ++jlane) {
                value(jblk, jlane) = next_value;
                next_value += Value{1};
            }
        }
    }
}

functionspace::BlockStructuredColumns make_block_fs(const StructuredGrid& grid, idx_t nproma) {
    util::Config config;
    config.set("halo", 0);
    config.set("nproma", nproma);
    return functionspace::BlockStructuredColumns(grid, config);
}

template <typename Value>
void check_relayout_field(const functionspace::BlockStructuredColumns& source_fs,
                          const functionspace::BlockStructuredColumns& target_fs,
                          const util::Config& options,
                          Value first_value) {
    Field source = source_fs.createField<Value>(option::name("source") | options);
    Field target = target_fs.createField<Value>(option::name("target") | options);

    fill_blocked_field(source_fs, source, first_value);
    prepare_source(source);
    prepare_target(target);

    copy_blocked_to_blocked(source, target, on_device());
    sync_target(target);

    EXPECT_EQ(source_fs.checksum(source), target_fs.checksum(target));
}

template <typename Value>
void check_blocked_to_nonblocked_field(const functionspace::StructuredColumns& structured_fs,
                                       const functionspace::BlockStructuredColumns& block_fs,
                                       const util::Config& options,
                                       Value first_value) {
    Field blocked = block_fs.createField<Value>(option::name("blocked") | options);
    Field target  = structured_fs.createField<Value>(option::name("target") | options);

    fill_blocked_field(block_fs, blocked, first_value);
    prepare_source(blocked);
    prepare_target(target);

    copy_blocked_to_nonblocked(blocked, target, on_device());
    sync_target(target);

    EXPECT_EQ(block_fs.checksum(blocked), structured_fs.checksum(target));
}

template <typename Value>
void check_nonblocked_to_blocked_field(const functionspace::StructuredColumns& structured_fs,
                                       const functionspace::BlockStructuredColumns& block_fs,
                                       const util::Config& options,
                                       Value first_value) {
    Field source = structured_fs.createField<Value>(option::name("source") | options);
    Field target = block_fs.createField<Value>(option::name("target") | options);

    fill_structured_field(source, first_value);
    prepare_source(source);
    prepare_target(target);

    copy_nonblocked_to_blocked(source, target, on_device());
    sync_target(target);

    EXPECT_EQ(structured_fs.checksum(source), block_fs.checksum(target));
}

}  // namespace

//-----------------------------------------------------------------------------

CASE("relayout copies blocked fields between different nproma layouts") {
    using Value = double;
    mpi::Scope scope("self");

    const auto grid       = StructuredGrid("O8");
    auto source_block_fs  = make_block_fs(grid, 5);
    auto target_block_fs  = make_block_fs(grid, 7);

    check_relayout_field<Value>(source_block_fs, target_block_fs, util::Config{}, Value{1});
    check_relayout_field<Value>(source_block_fs, target_block_fs, option::levels(4), Value{1000});
    check_relayout_field<Value>(source_block_fs, target_block_fs, option::variables(3) | option::levels(4), Value{2000});
}

CASE("relayout copies blocked fieldsets between different nproma layouts") {
    using Value = double;
    mpi::Scope scope("self");

    const auto grid       = StructuredGrid("O8");
    auto source_block_fs  = make_block_fs(grid, 5);
    auto target_block_fs  = make_block_fs(grid, 7);

    auto options_a = option::levels(3);
    auto options_b = option::variables(2) | option::levels(4);

    Field source_a = source_block_fs.createField<Value>(option::name("source_a") | options_a);
    Field source_b = source_block_fs.createField<Value>(option::name("source_b") | options_b);

    Field target_a = target_block_fs.createField<Value>(option::name("target_a") | options_a);
    Field target_b = target_block_fs.createField<Value>(option::name("target_b") | options_b);

    fill_blocked_field(source_block_fs, source_a, Value{3000});
    fill_blocked_field(source_block_fs, source_b, Value{4000});
    prepare_source(source_a);
    prepare_source(source_b);

    FieldSet source;
    source.add(source_a);
    source.add(source_b);

    FieldSet target;
    target.add(target_a);
    target.add(target_b);
    prepare_target(target);

    copy_blocked_to_blocked(source, target, on_device());
    sync_target(target);

    EXPECT_EQ(source_block_fs.checksum(source), target_block_fs.checksum(target));
}

CASE("relayout copies blocked fields to nonblocked fields") {
    using Value = double;
    mpi::Scope scope("self");

    const auto grid      = StructuredGrid("O8");
    auto structured_fs   = functionspace::StructuredColumns(grid, util::Config("halo", 0));
    auto block_fs        = make_block_fs(grid, 5);

    check_blocked_to_nonblocked_field<Value>(structured_fs, block_fs, util::Config{}, Value{5000});
    check_blocked_to_nonblocked_field<Value>(structured_fs, block_fs, option::levels(4), Value{6000});
    check_blocked_to_nonblocked_field<Value>(structured_fs, block_fs, option::variables(3) | option::levels(4), Value{7000});
}

CASE("relayout copies nonblocked fields to blocked fields") {
    using Value = double;
    mpi::Scope scope("self");

    const auto grid      = StructuredGrid("O8");
    auto structured_fs   = functionspace::StructuredColumns(grid, util::Config("halo", 0));
    auto block_fs        = make_block_fs(grid, 5);

    check_nonblocked_to_blocked_field<Value>(structured_fs, block_fs, util::Config{}, Value{8000});
    check_nonblocked_to_blocked_field<Value>(structured_fs, block_fs, option::levels(4), Value{9000});
    check_nonblocked_to_blocked_field<Value>(structured_fs, block_fs, option::variables(3) | option::levels(4), Value{10000});
}

CASE("relayout copies fieldsets between blocked and nonblocked layouts") {
    using Value = double;
    mpi::Scope scope("self");

    const auto grid      = StructuredGrid("O8");
    auto structured_fs   = functionspace::StructuredColumns(grid, util::Config("halo", 0));
    auto block_fs        = make_block_fs(grid, 5);

    auto options_a = option::levels(3);
    auto options_b = option::variables(2) | option::levels(4);

    Field nonblocked_source_a = structured_fs.createField<Value>(option::name("nonblocked_source_a") | options_a);
    Field nonblocked_source_b = structured_fs.createField<Value>(option::name("nonblocked_source_b") | options_b);
    EXPECT_EQ(nonblocked_source_a.rank(), 2);
    EXPECT_EQ(nonblocked_source_b.rank(), 3);

    Field blocked_source_a = block_fs.createField<Value>(option::name("blocked_source_a") | options_a);
    Field blocked_source_b = block_fs.createField<Value>(option::name("blocked_source_b") | options_b);
    EXPECT_EQ(blocked_source_a.rank(), 3);
    EXPECT_EQ(blocked_source_b.rank(), 4);
    
    Field nonblocked_target_a = structured_fs.createField<Value>(option::name("nonblocked_target_a") | options_a);
    Field nonblocked_target_b = structured_fs.createField<Value>(option::name("nonblocked_target_b") | options_b);
    EXPECT_EQ(nonblocked_target_a.rank(), 2);
    EXPECT_EQ(nonblocked_target_b.rank(), 3);

    Field blocked_target_a = block_fs.createField<Value>(option::name("blocked_target_a") | options_a);
    Field blocked_target_b = block_fs.createField<Value>(option::name("blocked_target_b") | options_b);
    EXPECT_EQ(blocked_target_a.rank(), 3);
    EXPECT_EQ(blocked_target_b.rank(), 4);

    fill_structured_field(nonblocked_source_a, Value{11000});
    fill_structured_field(nonblocked_source_b, Value{12000});
    prepare_source(nonblocked_source_a);
    prepare_source(nonblocked_source_b);

    fill_blocked_field(block_fs, blocked_source_a, Value{13000});
    fill_blocked_field(block_fs, blocked_source_b, Value{14000});
    prepare_source(blocked_source_a);
    prepare_source(blocked_source_b);

    FieldSet nonblocked_source;
    nonblocked_source.add(nonblocked_source_a);
    nonblocked_source.add(nonblocked_source_b);

    FieldSet blocked_source;
    blocked_source.add(blocked_source_a);
    blocked_source.add(blocked_source_b);

    FieldSet nonblocked_target;
    nonblocked_target.add(nonblocked_target_a);
    nonblocked_target.add(nonblocked_target_b);
    prepare_target(nonblocked_target);

    FieldSet blocked_target;
    blocked_target.add(blocked_target_a);
    blocked_target.add(blocked_target_b);
    prepare_target(blocked_target);

    copy_blocked_to_nonblocked(blocked_source, nonblocked_target, on_device());
    copy_nonblocked_to_blocked(nonblocked_source, blocked_target, on_device());
    sync_target(nonblocked_target);
    sync_target(blocked_target);

    EXPECT_EQ(block_fs.checksum(blocked_source), structured_fs.checksum(nonblocked_target));
    EXPECT_EQ(structured_fs.checksum(nonblocked_source), block_fs.checksum(blocked_target));
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
