/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/functionspace/BlockStructuredColumns.h"

#include <fstream>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/utils/MD5.h"

#include "atlas/array/Array.h"
#include "atlas/array/MakeView.h"
#include "atlas/domain.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/StructuredPartitionPolygon.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/fill.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Checksum.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/detail/Cache.h"

namespace atlas {
namespace functionspace {
namespace detail {

namespace {

template <class ValueType>
void block_copy(const Field sloc, Field loc, const functionspace::detail::BlockStructuredColumns& fs) {
    if (sloc.variables() and sloc.levels()) {
        auto loc_v  = array::make_view<ValueType, 4>(loc);
        auto sloc_v = array::make_view<ValueType, 3>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(2); ++jvar) {
                for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                    for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                        loc_v(jblk, jvar, jlev, jrof) = sloc_v(blk_begin+jrof, jlev, jvar);
                    }
                }
            }
        }
    }
    else if (not sloc.variables() and sloc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                    loc_v(jblk, jlev, jrof) = sloc_v(blk_begin+jrof, jlev);
                }
            }
        }
    }
    else if (sloc.variables() and not sloc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(1); ++jvar) {
                for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                    loc_v(jblk, jvar, jrof) = sloc_v(blk_begin+jrof, jvar);
                }
            }
        }
    }
    else {
        auto loc_v  = array::make_view<ValueType, 2>(loc);
        auto sloc_v = array::make_view<ValueType, 1>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                loc_v(jblk, jrof) = sloc_v(blk_begin+jrof);
            }
        }
    }
}

template <class ValueType>
void rev_block_copy(const Field loc, Field sloc, const functionspace::detail::BlockStructuredColumns& fs) {
    if (loc.variables() and loc.levels()) {
        auto loc_v  = array::make_view<ValueType, 4>(loc);
        auto sloc_v = array::make_view<ValueType, 3>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(2); ++jvar) {
                for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                    for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                        sloc_v(blk_begin+jrof, jlev, jvar) = loc_v(jblk, jvar, jlev, jrof);
                    }
                }
            }
        }
    }
    else if (not loc.variables() and loc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                    sloc_v(blk_begin+jrof, jlev) = loc_v(jblk, jlev, jrof);
                }
            }
        }
    }
    else if (loc.variables() and not loc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(1); ++jvar) {
                for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                    sloc_v(blk_begin+jrof, jvar) = loc_v(jblk, jvar, jrof);
                }
            }
        }
    }
    else {
        auto loc_v  = array::make_view<ValueType, 2>(loc);
        auto sloc_v = array::make_view<ValueType, 1>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            const idx_t blk_size = fs.block_size(jblk);
            const idx_t blk_begin = fs.block_begin(jblk);
            for (idx_t jrof = 0; jrof < blk_size; ++jrof) {
                sloc_v(blk_begin+jrof) = loc_v(jblk, jrof);
            }
        }
    }
}

void transpose_nonblocked_to_blocked(const Field& nonblocked, Field& blocked, const functionspace::detail::BlockStructuredColumns& fs) {
    auto kind = nonblocked.datatype().kind();
    if (kind == array::DataType::kind<int>()) {
        block_copy<int>(nonblocked, blocked, fs);
    }
    else if (kind == array::DataType::kind<long>()) {
        block_copy<long>(nonblocked, blocked, fs);
    }
    else if (kind == array::DataType::kind<unsigned long>()) {
        block_copy<unsigned long>(nonblocked, blocked, fs);
    }
    else if (kind == array::DataType::kind<float>()) {
        block_copy<float>(nonblocked, blocked, fs);
    }
    else if (kind == array::DataType::kind<double>()) {
        block_copy<double>(nonblocked, blocked, fs);
    }
    else {
        throw_Exception("datatype not supported", Here());
    }
}

void transpose_blocked_to_nonblocked(const Field& blocked, Field& nonblocked, const functionspace::detail::BlockStructuredColumns& fs) {
    auto kind = blocked.datatype().kind();
    if (kind == array::DataType::kind<int>()) {
        rev_block_copy<int>(blocked, nonblocked, fs);
    }
    else if (kind == array::DataType::kind<long>()) {
        rev_block_copy<long>(blocked, nonblocked, fs);
    }
    else if (kind == array::DataType::kind<unsigned long>()) {
        rev_block_copy<unsigned long>(blocked, nonblocked, fs);
    }
    else if (kind == array::DataType::kind<float>()) {
        rev_block_copy<float>(blocked, nonblocked, fs);
    }
    else if (kind == array::DataType::kind<double>()) {
        rev_block_copy<double>(blocked, nonblocked, fs);
    }
    else {
        throw_Exception("datatype not supported", Here());
    }
}


}// namespace

array::ArrayShape BlockStructuredColumns::config_shape(const eckit::Configuration& config) const {
    array::ArrayShape shape;

    bool global = false;
    config.get("global", global);
    if (global) {
        return structuredcolumns_->config_shape(config);
    }
    else {
        shape.emplace_back(nblks_);
        idx_t variables(0);
        config.get("variables", variables);
        if (variables > 0) {
            shape.emplace_back(variables);
        }
        idx_t levels(structuredcolumns_->levels());
        config.get("levels", levels);
        if (levels > 0) {
            shape.emplace_back(levels);
        }
        shape.emplace_back(nproma_);
    }
    return shape;
}

array::ArraySpec BlockStructuredColumns::config_spec(const eckit::Configuration& config) const {
    return array::ArraySpec(config_shape(config), structuredcolumns_->config_alignment(config));
}

// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const eckit::Configuration& config):
    BlockStructuredColumns::BlockStructuredColumns(grid, grid::Partitioner(), config) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const grid::Partitioner& p, const eckit::Configuration& config):
    BlockStructuredColumns(grid, Vertical(config), p, config) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const grid::Distribution& distribution,
                                     const eckit::Configuration& config):
    BlockStructuredColumns(grid, distribution, Vertical(config), config) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const grid::Distribution& distribution, const Vertical& vertical,
                                     const eckit::Configuration& config):
    structuredcolumns_(new StructuredColumns(grid, distribution, vertical, config)),
    structuredcolumns_handle_(structuredcolumns_){
    setup(config);
}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const Vertical& vertical, const eckit::Configuration& config):
    BlockStructuredColumns(grid, vertical, grid::Partitioner(), config) {}

BlockStructuredColumns::BlockStructuredColumns(const Grid& grid, const Vertical& vertical, const grid::Partitioner& p,
                                     const eckit::Configuration& config):
    structuredcolumns_(new StructuredColumns(grid, vertical, p, config)),
    structuredcolumns_handle_(structuredcolumns_){
    setup(config);
}

// ----------------------------------------------------------------------------
// Create Field
// ----------------------------------------------------------------------------
Field BlockStructuredColumns::createField(const eckit::Configuration& options) const {
    Field field(structuredcolumns_->config_name(options), structuredcolumns_->config_datatype(options), config_spec(options));
    structuredcolumns_->set_field_metadata(options, field);
    field.set_functionspace(this);

    bool global = false;
    options.get("global", global);
    if (not global) {
        field.set_horizontal_dimension({0,field.rank()-1});
    }

    return field;
}

Field BlockStructuredColumns::createField(const Field& other, const eckit::Configuration& config) const {
    return createField(option::name(other.name()) | option::datatype(other.datatype()) |
                       option::levels(other.levels()) |option::variables(other.variables()) |
                       option::type(other.metadata().getString("type", "scalar")) | config);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Scatter FieldSet
// ----------------------------------------------------------------------------
void BlockStructuredColumns::scatter(const FieldSet& global_fieldset, FieldSet& local_fieldset) const {
    ATLAS_ASSERT(local_fieldset.size() == global_fieldset.size());
    for (idx_t f = 0; f < local_fieldset.size(); ++f) {
        const Field& glb  = global_fieldset[f];
        Field& loc        = local_fieldset[f];
        auto sloc         = structuredcolumns_->createField(glb, util::Config("global",false));
        structuredcolumns_->scatter(glb, sloc);
        loc.metadata() = sloc.metadata();
        transpose_nonblocked_to_blocked(sloc, loc, *this);
    }
}

// ----------------------------------------------------------------------------
// Scatter Field
// ----------------------------------------------------------------------------
void BlockStructuredColumns::scatter(const Field& global, Field& local) const {
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add(global);
    local_fields.add(local);
    scatter(global_fields, local_fields);
}

// ----------------------------------------------------------------------------
// Gather FieldSet
// ----------------------------------------------------------------------------
void BlockStructuredColumns::gather(const FieldSet& local_fieldset, FieldSet& global_fieldset) const {
    ATLAS_ASSERT(local_fieldset.size() == global_fieldset.size());
    for (idx_t f = 0; f < local_fieldset.size(); ++f) {
        const Field& loc = local_fieldset[f];
        Field& glb       = global_fieldset[f];
        auto sloc = structuredcolumns_->createField(loc, util::Config("global",false));
        transpose_blocked_to_nonblocked(loc, sloc, *this);
        structuredcolumns_->gather(sloc, glb);
        glb.metadata() = loc.metadata();
        glb.metadata().set("global", false);
    }
}
// ----------------------------------------------------------------------------

void BlockStructuredColumns::setup(const eckit::Configuration &config) {
    nproma_ = 1;
    idx_t tmp_nproma;
    if (config.get("nproma", tmp_nproma)) {
        ATLAS_ASSERT(tmp_nproma > 0);
        nproma_ = tmp_nproma;
    }

    nblks_ = std::floor( structuredcolumns_->size() / nproma_ );
    endblk_size_ = nproma_;
    if (structuredcolumns_->size() % nproma_ > 0) {
        endblk_size_ = structuredcolumns_->size() - nblks_ * nproma_;
        nblks_++;
    }
}

// ----------------------------------------------------------------------------
// Gather Field
// ----------------------------------------------------------------------------
void BlockStructuredColumns::gather(const Field& local, Field& global) const {
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add(local);
    global_fields.add(global);
    gather(local_fields, global_fields);
}

// ----------------------------------------------------------------------------
// Checksum Field
// ----------------------------------------------------------------------------

namespace {

struct ChecksumData {
    Field local_lane_checksums;
    Field global_lane_checksums;
    int root{0};
};

bool block_slice_contiguous(const Field& f) {
    idx_t expected_stride = 1;
    for (idx_t dim = f.rank() - 1; dim > 0; --dim) {
        if (f.stride(dim) != expected_stride) {
            return false;
        }
        expected_stride *= f.shape(dim);
    }
    return true;
}

template <typename T>
void checksum_update_lane(util::checksum_t& checksum, const array::Array& f, const T* data, idx_t dim, idx_t offset) {
    if (dim == f.rank() - 1) {
        util::checksum_update(checksum, data + offset, 1);
        return;
    }
    for (idx_t index = 0; index < f.shape(dim); ++index) {
        checksum_update_lane(checksum, f, data, dim + 1, offset + index * f.stride(dim));
    }
}

template <typename T>
void update_lane_checksums_with_field_contiguous(const BlockStructuredColumns& fs, const Field& f, const T* data,
                                                 array::ArrayView<util::checksum_t, 2>& local_lane_checksums_view) {
    const idx_t nproma = f.shape(f.rank() - 1);
    const idx_t full_blocks = (fs.nblks() > 0 ? fs.nblks() - 1 : 0);

    if (f.rank() == 4) {
        auto update_block = [&](idx_t jblk, idx_t block_size) {
            const idx_t block_offset = jblk * f.stride(0);
            const T* block_data = data + block_offset;
            // [nblks, nvar, nlev, nproma]: iterate lev-outer, var-inner
            // to match global field order [npoints, nlev, nvar]
            const idx_t nvar = f.shape(1);
            const idx_t nlev = f.shape(2);
            const size_t lane_size = static_cast<size_t>(nvar);
            const size_t lane_stride = static_cast<size_t>(nlev * nproma);
            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                for (idx_t jlane = 0; jlane < block_size; ++jlane) {
                    const T* lane_data = block_data + jlev * nproma + jlane;
                    util::checksum_update(local_lane_checksums_view(jblk, jlane), lane_data, lane_size, lane_stride);
                }
            }
        };

        atlas_omp_parallel_for(idx_t jblk = 0; jblk < full_blocks; ++jblk) {
            update_block(jblk, nproma);
        }
        if (fs.nblks() > 0) {
            update_block(fs.nblks() - 1, fs.block_size(fs.nblks() - 1));
        }
    }
    else {
        const size_t lane_stride = static_cast<size_t>(nproma);
        size_t lane_size{1};
        for (idx_t dim = 1; dim < f.rank() - 1; ++dim) {
            lane_size *= f.shape(dim);
        }
        auto update_block = [&](idx_t jblk, idx_t block_size) {
            const idx_t block_offset = jblk * f.stride(0);
            const T* block_data = data + block_offset;
            for (idx_t jlane = 0; jlane < block_size; ++jlane) {
                const T* lane_data = block_data + jlane;
                util::checksum_update(local_lane_checksums_view(jblk, jlane), lane_data, lane_size, lane_stride);
            }
        };

        atlas_omp_parallel_for(idx_t jblk = 0; jblk < full_blocks; ++jblk) {
            update_block(jblk, nproma);
        }
        if (fs.nblks() > 0) {
            update_block(fs.nblks() - 1, fs.block_size(fs.nblks() - 1));
        }
    }
}

template <typename T>
void update_lane_checksums_with_field(const BlockStructuredColumns& fs, const Field& f,
                                      ChecksumData& checksum_data) {
    auto local_lane_checksums_view = array::make_view<util::checksum_t, 2>(checksum_data.local_lane_checksums);
    const auto* data = f.array().host_data<T>();
    if (block_slice_contiguous(f)) {
        update_lane_checksums_with_field_contiguous(fs, f, data, local_lane_checksums_view);
        return;
    }

    // A fallback; we should usually not be here normally, as a field block is usually contiguous
    const idx_t lane_stride = f.stride(f.rank() - 1);
    atlas_omp_parallel_for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
        const idx_t block_offset = jblk * f.stride(0);
        const idx_t block_size   = fs.block_size(jblk);
        if (f.rank() == 4) {
            // Iterate lev-outer, var-inner to match global field order [npoints, nlev, nvar]
            const idx_t nvar = f.shape(1);
            const idx_t nlev = f.shape(2);
            for (idx_t jlane = 0; jlane < block_size; ++jlane) {
                util::checksum_t lane_checksum = local_lane_checksums_view(jblk, jlane);
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                        const idx_t offset = block_offset + jlane * lane_stride + jvar * f.stride(1) + jlev * f.stride(2);
                        util::checksum_update(lane_checksum, data + offset, size_t{1});
                    }
                }
                local_lane_checksums_view(jblk, jlane) = lane_checksum;
            }
        }
        else {
            for (idx_t jlane = 0; jlane < block_size; ++jlane) {
                util::checksum_t lane_checksum = local_lane_checksums_view(jblk, jlane);
                checksum_update_lane(lane_checksum, f, data, idx_t{1}, block_offset + jlane * lane_stride);
                local_lane_checksums_view(jblk, jlane) = lane_checksum;
            }
        }
    }
}

void update_lane_checksums_with_field_dynamic(const BlockStructuredColumns& fs, const Field& f,
                                              ChecksumData& checksum_data) {
    switch (f.datatype().kind()) {
        case array::DataType::kind<int>():    return update_lane_checksums_with_field<int>(fs, f, checksum_data);
        case array::DataType::kind<long>():   return update_lane_checksums_with_field<long>(fs, f, checksum_data);
        case array::DataType::kind<float>():  return update_lane_checksums_with_field<float>(fs, f, checksum_data);
        case array::DataType::kind<double>(): return update_lane_checksums_with_field<double>(fs, f, checksum_data);
        default:
            throw_Exception("datatype not supported", Here());
    }
}

ChecksumData get_checksum_data(const BlockStructuredColumns& fs) {
    pluto::scope scope;
    pluto::host::set_default_resource(pluto::host_pool_resource());
    ChecksumData checksum_data;
    checksum_data.root = 0;
    checksum_data.local_lane_checksums = fs.createField(option::name("checksum_local_lane_checksums") |
                                                        option::datatypeT<util::checksum_t>() |
                                                        option::levels(0));
    checksum_data.global_lane_checksums = fs.createField(checksum_data.local_lane_checksums, option::global(checksum_data.root));
    auto local_lane_checksums_view = array::make_view<util::checksum_t, 2>(checksum_data.local_lane_checksums);
    local_lane_checksums_view.assign(0);
    return checksum_data;
}

util::checksum_t checksum_from_lane_states(const BlockStructuredColumns& fs, ChecksumData& checksum_data) {
    auto& local_lane_checksums = checksum_data.local_lane_checksums;
    auto local_lane_checksums_view = array::make_view<util::checksum_t, 2>(local_lane_checksums);
    for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
        for (idx_t jlane = 0; jlane < fs.nproma(); ++jlane) {
            local_lane_checksums_view(jblk, jlane) = util::checksum_digest(local_lane_checksums_view(jblk, jlane));
        }
    }

    auto& mpi_comm = mpi::comm(fs.mpi_comm());

    if (mpi_comm.size() == 1) {
        ATLAS_ASSERT(local_lane_checksums.contiguous());
        const auto* local_lane_checksums_data = local_lane_checksums.array().host_data<util::checksum_t>();
        return util::checksum(local_lane_checksums_data, static_cast<size_t>(fs.structuredcolumns().sizeOwned()));
    }

    auto& global_lane_checksums = checksum_data.global_lane_checksums;
    fs.gather(local_lane_checksums, global_lane_checksums);

    util::checksum_t global_checksum = 0;
    if (mpi_comm.rank() == checksum_data.root) {
        ATLAS_ASSERT(global_lane_checksums.contiguous());
        const auto* global_lane_checksums_data = global_lane_checksums.array().host_data<util::checksum_t>();
        global_checksum = util::checksum(global_lane_checksums_data, global_lane_checksums.size());
    }
    mpi_comm.broadcast(global_checksum, checksum_data.root);
    return global_checksum;
}

util::checksum_t field_checksum(const BlockStructuredColumns& fs, const Field& f) {
    ATLAS_ASSERT(f.rank() > 1);
    auto checksum_data = get_checksum_data(fs);
    update_lane_checksums_with_field_dynamic(fs, f, checksum_data);
    return checksum_from_lane_states(fs, checksum_data);
}

util::checksum_t fieldset_checksum(const BlockStructuredColumns& fs, const FieldSet& fields) {
    auto checksum_data = get_checksum_data(fs);
    for (idx_t jfld = 0; jfld < fields.size(); ++jfld) {
        update_lane_checksums_with_field_dynamic(fs, fields[jfld], checksum_data);
    }
    return checksum_from_lane_states(fs, checksum_data);
}
}  // namespace


std::string BlockStructuredColumns::checksum(const Field& f) const {
    ATLAS_TRACE("BlockStructuredColumns::checksum");
    return util::checksum_to_hex_str(field_checksum(*this, f));
}

std::string BlockStructuredColumns::checksum(const FieldSet& f) const {
    ATLAS_TRACE("BlockStructuredColumns::checksum");
    return util::checksum_to_hex_str(fieldset_checksum(*this, f));
}

// ----------------------------------------------------------------------------


}  // namespace detail
}  // namespace functionspace
}  // namespace atlas
