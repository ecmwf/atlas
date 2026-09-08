/*
 * (C) Copyright 2026- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

#include "pluto/pluto.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Checksum.h"

namespace atlas {
namespace functionspace {
namespace detail {

/// @brief Shared checksum helper for column-based function spaces.
///
/// This helper implements the lane-based checksum algorithm used by
/// NodeColumns, EdgeColumns, and CellColumns. It is an internal utility in the
/// functionspace detail layer and is not part of the public Atlas API.
///
/// The class assumes the function space exposes the minimal interface needed by
/// the algorithm: `size()`, `ghost()`, `gather()`, and `mpi_comm()`.
///
/// @tparam FunctionSpace Internal function space implementation type.
template <typename FunctionSpace>
class ColumnChecksum {
public:
    /// @brief Construct a checksum helper for one function space instance.
    ///
    /// @param functionspace Internal function space instance providing entity
    ///                      ownership and gather operations.
    /// @param global_size   Global number of entities represented on the root
    ///                      rank after gather.
    ColumnChecksum(const FunctionSpace& functionspace, size_t global_size): functionspace_(functionspace), global_size_(global_size) {}

    /// @brief Compute the checksum of a field set.
    ///
    /// The checksum is built from owned entities only, then gathered and
    /// reduced to a global hexadecimal digest.
    ///
    /// @param fieldset Field set defined on the associated function space.
    /// @return Hexadecimal checksum string.
    std::string compute(const FieldSet& fieldset) const {
        return util::checksum_to_hex_str(compute_checksum(fieldset));
    }

    /// @brief Compute the checksum of a single field.
    ///
    /// This is a convenience overload that wraps the field into a temporary
    /// FieldSet and reuses the FieldSet-based implementation.
    ///
    /// @param field Field defined on the associated function space.
    /// @return Hexadecimal checksum string made of 4 characters (16 bits), see util::checksum_to_hex_str.
    std::string compute(const Field& field) const {
        FieldSet fieldset;
        fieldset.add(field);
        return compute(fieldset);
    }

private:
    template <typename Value>
    static void checksum_update_point(util::checksum_t& checksum, const array::Array& field, const Value* data, idx_t dim,
                                      idx_t offset) {
        if (dim == field.rank()) {
            util::checksum_update(checksum, data + offset, size_t{1});
            return;
        }
        for (idx_t index = 0; index < field.shape(dim); ++index) {
            checksum_update_point(checksum, field, data, dim + 1, offset + index * field.stride(dim));
        }
    }

    static bool is_owned(const array::ArrayView<int, 1>& ghost, idx_t index) {
        return ghost(index) == 0;
    }

    static idx_t owned_contiguous_prefix_size(const array::ArrayView<int, 1>& ghost, idx_t size) {
        idx_t owned_size = size;
        while (owned_size > 0 && ghost(owned_size - 1) != 0) {
            --owned_size;
        }
        for (idx_t index = 0; index < owned_size; ++index) {
            if (ghost(index) != 0) {
                return -1;
            }
        }
        return owned_size;
    }

    template <typename Value>
    void update_lane_checksums_with_field(const Field& field, const array::ArrayView<int, 1>& ghost,
                                          util::checksum_t lane_checksums[], size_t lane_checksums_size) const {
        const idx_t npts = std::min(functionspace_.size(), field.shape(0));
        ATLAS_ASSERT(lane_checksums_size >= static_cast<size_t>(npts));

        const auto* data = field.array().host_data<Value>();

        if (field.contiguous()) {
            const auto point_stride     = field.stride(0);
            const idx_t contiguous_npts = owned_contiguous_prefix_size(ghost, functionspace_.size());
            if (contiguous_npts >= 0) {
                atlas_omp_parallel_for(idx_t index = 0; index < contiguous_npts; ++index) {
                    util::checksum_update(lane_checksums[index], data + index * point_stride, point_stride);
                }
            }
            else {
                atlas_omp_parallel_for(idx_t index = 0; index < npts; ++index) {
                    if (is_owned(ghost, index)) {
                        util::checksum_update(lane_checksums[index], data + index * point_stride, point_stride);
                    }
                }
            }
            return;
        }

        atlas_omp_parallel_for(idx_t index = 0; index < npts; ++index) {
            if (is_owned(ghost, index)) {
                util::checksum_t lane = lane_checksums[index];
                checksum_update_point(lane, field, data, idx_t{1}, index * field.stride(0));
                lane_checksums[index] = lane;
            }
        }
    }

    void update_lane_checksums_with_field_dynamic(const Field& field, const array::ArrayView<int, 1>& ghost,
                                                  util::checksum_t lane_checksums[], size_t lane_checksums_size) const {
        switch (field.datatype().kind()) {
            case array::DataType::kind<int>():
                return update_lane_checksums_with_field<int>(field, ghost, lane_checksums, lane_checksums_size);
            case array::DataType::kind<long>():
                return update_lane_checksums_with_field<long>(field, ghost, lane_checksums, lane_checksums_size);
            case array::DataType::kind<float>():
                return update_lane_checksums_with_field<float>(field, ghost, lane_checksums, lane_checksums_size);
            case array::DataType::kind<double>():
                return update_lane_checksums_with_field<double>(field, ghost, lane_checksums, lane_checksums_size);
            default:
                throw_Exception("datatype not supported", Here());
        }
    }

    util::checksum_t compute_checksum(const FieldSet& fieldset) const {
        pluto::scope scope;
        auto* memory_resource = pluto::host_pool_resource();
        pluto::host::set_default_resource(memory_resource);

        const idx_t npts = functionspace_.size();
        Field ghost_field = functionspace_.ghost();
        auto ghost        = array::make_view<int, 1>(ghost_field);
        std::vector<util::checksum_t, pluto::allocator<util::checksum_t>> local_lanes(npts, 0, memory_resource);
        for (idx_t field_idx = 0; field_idx < fieldset.size(); ++field_idx) {
            update_lane_checksums_with_field_dynamic(fieldset[field_idx], ghost, local_lanes.data(), local_lanes.size());
        }

        for (idx_t index = 0; index < npts; ++index) {
            if (is_owned(ghost, index)) {
                local_lanes[index] = util::checksum_digest(local_lanes[index]);
            }
        }

        const idx_t root = 0;
        auto& comm       = mpi::comm(functionspace_.mpi_comm());
        size_t nb_global = (comm.rank() == root) ? global_size_ : 0;

        std::unique_ptr<util::checksum_t[], std::function<void(util::checksum_t*)>> global_lanes(
            nb_global > 0 ? pluto::allocator<util::checksum_t>(memory_resource).allocate(nb_global) : nullptr,
            [&](util::checksum_t* ptr) {
                if (ptr) {
                    pluto::allocator<util::checksum_t>(memory_resource).deallocate(ptr, nb_global);
                }
            });

        {
            mpi::Scope mpi_scope(functionspace_.mpi_comm());
            idx_t lvar_stride = 1;
            idx_t lvar_shape  = 1;
            idx_t gvar_stride = 1;
            idx_t gvar_shape  = 1;
            functionspace_.gather().gather(local_lanes.data(), &lvar_stride, &lvar_shape, 1, global_lanes.get(),
                                           &gvar_stride, &gvar_shape, 1, root);
        }

        util::checksum_t global_checksum = 0;
        if (comm.rank() == root) {
            global_checksum = util::checksum(global_lanes.get(), nb_global);
        }
        comm.broadcast(global_checksum, root);
        return global_checksum;
    }

private:
    const FunctionSpace& functionspace_;
    size_t global_size_;
};

}  // namespace detail
}  // namespace functionspace
}  // namespace atlas