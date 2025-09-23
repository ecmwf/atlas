/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "AssembleGlobalMatrix.h"

#include <unordered_map>

#include "eckit/linalg/types.h"

#include "atlas/array.h"
#include "atlas/linalg/sparse/SparseMatrixToTriplets.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/functionspace/NodeColumns.h"

namespace atlas::interpolation {

linalg::SparseMatrixStorage assemble_global_matrix(const Interpolation& interpolation, int mpi_root) {

    auto src_fs = interpolation.source();
    auto tgt_fs = interpolation.target();

    auto& mpi_comm = mpi::comm();
    auto mpi_size  = mpi_comm.size();
    auto mpi_rank  = mpi_comm.rank();

    std::vector<gidx_t> global_cols;
    std::vector<gidx_t> global_rows;
    std::vector<double> global_vals;

    // Compute global_cols, global_rows, global_vals
    {
        std::vector<gidx_t> local_cols;
        std::vector<gidx_t> local_rows;
        std::vector<double> local_vals;

        // Compute local_cols, local_rows, local_vals
        {
            interpolation::MatrixCache interpolation_cache = interpolation.createCache();
            auto local_matrix_storage = interpolation_cache.matrix();
            auto local_matrix = atlas::linalg::make_host_view<eckit::linalg::Scalar,eckit::linalg::Index>(local_matrix_storage);

            const auto src_ridx = array::make_indexview<idx_t, 1>(src_fs.remote_index());

            const auto src_global_index = array::make_view<gidx_t, 1>(src_fs.global_index());
            const auto tgt_global_index = array::make_view<gidx_t, 1>(tgt_fs.global_index());
            const auto src_part         = array::make_view<int, 1>(src_fs.partition());
            const auto tgt_ghost        = array::make_view<int, 1>(tgt_fs.ghost());

            eckit::mpi::Buffer<gidx_t> recv_global_idx_buf(mpi_size);
            {
                std::vector<gidx_t> send_global_idx(src_global_index.size());
                for (int i = 0; i < send_global_idx.size(); ++i) {
                    send_global_idx[i] = src_global_index[i];
                }
                mpi_comm.allGatherv(send_global_idx.begin(), send_global_idx.end(), recv_global_idx_buf);
            }
            auto recv_global_idx = [&recv_global_idx_buf](int part, int ridx) {
                return recv_global_idx_buf.buffer[ recv_global_idx_buf.displs[part] + ridx ];
            };
            auto src_gidx = [&](auto idx) {
                return recv_global_idx(src_part(idx), src_ridx(idx));
            };

            local_cols.reserve(local_matrix.nnz());
            local_rows.reserve(local_matrix.nnz());
            local_vals.reserve(local_matrix.nnz());

            const auto* outer  = local_matrix.outer();
            const auto* inner  = local_matrix.inner();
            const auto* value  = local_matrix.value();
            const size_t nrows = local_matrix.rows();

            for(idx_t r = 0; r < nrows; ++r) {
                if( not tgt_ghost(r) ) {
                    idx_t row = tgt_global_index(r);
                    for (idx_t c = outer[r]; c < outer[r + 1]; ++c) {
                        auto src_idx = inner[c];
                        auto col = src_gidx(src_idx);
                        auto val = value[c];
                        local_rows.emplace_back(row);
                        local_cols.emplace_back(col);
                        local_vals.emplace_back(val);
                    }
                }
            }
        }

        // Gather local to global on mpi_root
        size_t local_nnz = local_vals.size();
        std::vector<size_t> local_nnz_per_rank(mpi_size);
        mpi_comm.gather(local_nnz, local_nnz_per_rank, mpi_root);

        constexpr int mpi_tag = 0;
        if (mpi_rank == mpi_root) {
            size_t global_nnz = std::accumulate(local_nnz_per_rank.begin(), local_nnz_per_rank.end(), 0);
            global_cols.resize(global_nnz);
            global_rows.resize(global_nnz);
            global_vals.resize(global_nnz);
            size_t pos = 0;
            size_t pos_root = 0;
            for (int jproc = 0; jproc < mpi_size; ++jproc) {
                if (jproc == mpi_root) {
                    pos += (jproc == 0) ? 0 : local_nnz_per_rank[jproc - 1];
                    pos_root = pos;
                    continue;
                }
                pos += (jproc == 0 ? 0 : local_nnz_per_rank[jproc - 1]);
                mpi_comm.receive(global_cols.data() + pos, local_nnz_per_rank[jproc], jproc, mpi_tag);
                mpi_comm.receive(global_rows.data() + pos, local_nnz_per_rank[jproc], jproc, mpi_tag);
                mpi_comm.receive(global_vals.data() + pos, local_nnz_per_rank[jproc], jproc, mpi_tag);
            }
            for (int i = 0; i < local_nnz; ++i) {
                global_cols[i + pos_root] = local_cols[i];
                global_rows[i + pos_root] = local_rows[i];
                global_vals[i + pos_root] = local_vals[i];
            }
        }
        else {
            mpi_comm.send(local_cols.data(), local_cols.size(), mpi_root, mpi_tag);
            mpi_comm.send(local_rows.data(), local_rows.size(), mpi_root, mpi_tag);
            mpi_comm.send(local_vals.data(), local_vals.size(), mpi_root, mpi_tag);
        }
    }


    auto compute_max_global_index = [&mpi_comm,mpi_root](const FunctionSpace& fs) {
        auto global_index = array::make_view<gidx_t, 1>(fs.global_index());
        auto ghost        = array::make_view<int, 1>(fs.ghost());
        gidx_t max_gidx{0};
        for( size_t i = 0; i < global_index.size(); ++i) {
            if (not ghost(i)) {
                max_gidx = std::max(max_gidx, global_index(i));
            }
        }
#if ATLAS_ECKIT_VERSION_AT_LEAST(1, 22, 2)
        // eckit::mpi::Comm::reduceInPlace() was only added in eckit 1.22.2
        mpi_comm.reduceInPlace(max_gidx, eckit::mpi::max(), mpi_root);
#else
        mpi_comm.allReduceInPlace(max_gidx, eckit::mpi::max());
        (void)mpi_root; // to silence unused lambda capture warning
#endif
        return max_gidx;
    };

    // compute max column and row indices
    gidx_t tgt_max_gidx = compute_max_global_index(tgt_fs);
    gidx_t src_max_gidx = compute_max_global_index(src_fs);

    linalg::SparseMatrixStorage global_matrix;
    if (mpi_rank == mpi_root) {
        size_t nrows = tgt_max_gidx;
        size_t ncols = src_max_gidx;
        size_t index_base = 1;
        bool is_sorted = false;
        global_matrix = atlas::linalg::make_sparse_matrix_storage_from_rows_columns_values(nrows, ncols, global_rows, global_cols, global_vals, index_base, is_sorted);
    }
    return global_matrix;
}

template <typename ViewValue, typename ViewIndex, typename Value, typename Index> 
void distribute_global_matrix(const linalg::SparseMatrixView<ViewValue,ViewIndex>& global_matrix,
    const array::Array& partition, std::vector<Index>& rows, std::vector<Index>& cols, std::vector<Value>& vals, int mpi_root) {
    ATLAS_TRACE("distribute_global_matrix_lowlevel");

    const auto tgt_part_glb = array::make_view<int,1>(partition);
    
    auto& mpi_comm = mpi::comm();
    auto mpi_size  = mpi_comm.size();
    auto mpi_rank  = mpi_comm.rank();

    // compute how many nnz-entries each task gets
    size_t nnz_loc = 0;
    int mpi_tag = 0;
    std::vector<std::size_t> nnz_per_task(mpi_size);
    const auto outer  = global_matrix.outer();
    const auto inner  = global_matrix.inner();
    const auto value  = global_matrix.value();
    if (mpi_rank == mpi_root) {

        ATLAS_ASSERT(global_matrix.rows() <= tgt_part_glb.size());
        for(std::size_t r = 0; r < global_matrix.rows(); ++r) {
            if (tgt_part_glb(r) >= mpi_size) {
                    ATLAS_DEBUG_VAR(r);
                    ATLAS_DEBUG_VAR(mpi_size);
                    ATLAS_DEBUG_VAR(tgt_part_glb(r));
                    ATLAS_ASSERT(tgt_part_glb(r) < nnz_per_task.size());
            }
            nnz_per_task[tgt_part_glb(r)] += outer[r+1] - outer[r];
        }
        for (int jproc = 0; jproc < mpi::comm().size(); ++jproc) {
            if (jproc != mpi_root) {
                mpi::comm().send(nnz_per_task.data() + jproc, 1, jproc, mpi_tag);
            }
        }
        nnz_loc = nnz_per_task[mpi_root];
    }
    else {
        mpi_comm.receive(&nnz_loc, 1, mpi_root, mpi_tag);
    }

    rows.resize(nnz_loc);
    cols.resize(nnz_loc);
    vals.resize(nnz_loc);

    if (mpi_rank == mpi_root) {
        std::vector<std::vector<Index>> send_rows(mpi_size);
        std::vector<std::vector<Index>> send_cols(mpi_size);
        std::vector<std::vector<Value>> send_vals(mpi_size);
        for(std::size_t jproc=0; jproc < mpi_size; ++jproc) {
            send_rows[jproc].reserve(nnz_per_task[jproc]);
            send_cols[jproc].reserve(nnz_per_task[jproc]);
            send_vals[jproc].reserve(nnz_per_task[jproc]);
        }
        for(std::size_t r = 0; r < global_matrix.rows(); ++r) {
            int jproc = tgt_part_glb(r);
            for (auto c = outer[r]; c < outer[r + 1]; ++c) {
                auto col = inner[c];
                send_rows[jproc].emplace_back(r);
                send_cols[jproc].emplace_back(col);
                send_vals[jproc].emplace_back(value[c]);
            }
        }
        for(std::size_t jproc = 0; jproc < mpi_size; ++jproc) {
            if (jproc != mpi_root) {
                mpi_comm.send(send_rows[jproc].data(), send_rows[jproc].size(), jproc, mpi_tag);
                mpi_comm.send(send_cols[jproc].data(), send_cols[jproc].size(), jproc, mpi_tag);
                mpi_comm.send(send_vals[jproc].data(), send_vals[jproc].size(), jproc, mpi_tag);
            }
            else {
                rows = send_rows[jproc];
                cols = send_cols[jproc];
                vals = send_vals[jproc];
            }
        }
    }
    else {
        mpi_comm.receive(rows.data(), nnz_loc, mpi_root, mpi_tag);
        mpi_comm.receive(cols.data(), nnz_loc, mpi_root, mpi_tag);
        mpi_comm.receive(vals.data(), nnz_loc, mpi_root, mpi_tag);
    }
}

linalg::SparseMatrixStorage distribute_global_matrix(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs, const linalg::SparseMatrixStorage& gmatrix, int mpi_root) {
    ATLAS_TRACE("distribute_global_matrix");
    Field field_tgt_part_glb = tgt_fs.createField(tgt_fs.partition(), option::global(mpi_root));

    // This should not be needed, but there's a problem with the NodeColumns functionspace for ORCA grids leading to uninitialized partitions
    if (functionspace::NodeColumns(tgt_fs)) {
        if (auto grid = functionspace::NodeColumns(tgt_fs).mesh().grid()) {
            if (grid.type() == "ORCA") {
                array::make_view<int,1>(field_tgt_part_glb).assign(0);
            }
        }
    }

    ATLAS_TRACE_SCOPE("gather partition") {
        tgt_fs.gather(tgt_fs.partition(), field_tgt_part_glb);
    }
    using Index = eckit::linalg::Index;
    using Value = eckit::linalg::Scalar;
    std::vector<Index> rows, cols;
    std::vector<Value> vals;
    distribute_global_matrix(atlas::linalg::make_host_view<Value, Index>(gmatrix), field_tgt_part_glb, rows, cols, vals, mpi_root);

    // map global index to local index
    std::unordered_map<gidx_t, idx_t> to_local_rows;
    std::unordered_map<gidx_t, idx_t> to_local_cols;

    ATLAS_TRACE_SCOPE("convert to local indexing") {
        auto tgt_gidx_exchanged = tgt_fs.createField(tgt_fs.global_index());
        tgt_gidx_exchanged.array().copy(tgt_fs.global_index());
        tgt_fs.haloExchange(tgt_gidx_exchanged);
        const auto tgt_global_index = array::make_view<gidx_t, 1>(tgt_gidx_exchanged);
        const auto tgt_ghost = array::make_view<int,1>(tgt_fs.ghost());

        auto src_gidx_exchanged = src_fs.createField(src_fs.global_index());
        src_gidx_exchanged.array().copy(src_fs.global_index());
        src_fs.haloExchange(src_gidx_exchanged);
        const auto src_global_index = array::make_view<gidx_t, 1>(src_gidx_exchanged);
        const auto src_ghost = array::make_view<int,1>(src_fs.ghost());

        for (idx_t r = 0; r < tgt_global_index.size(); ++r) {
            auto gr = tgt_global_index(r);
            if (tgt_ghost(r) && to_local_rows.find(gr) != to_local_rows.end()) {
                continue;
            }
            to_local_rows[gr] = r;
        }
        for (idx_t c = 0; c < src_global_index.size(); ++c) {
            auto gc = src_global_index(c);
            if (src_ghost(c) && to_local_cols.find(gc) != to_local_cols.end()) {
                continue;
            }
            to_local_cols[gc] = c;
        }
        for (int r = 0; r < rows.size(); ++r) {
            rows[r] = to_local_rows[rows[r] + 1];
            cols[r] = to_local_cols[cols[r] + 1];
        }
    }

    linalg::SparseMatrixStorage matrix;
    constexpr int index_base = 0;
    constexpr bool is_sorted = false;
    matrix = linalg::make_sparse_matrix_storage_from_rows_columns_values(tgt_fs.size(), src_fs.size(), rows, cols, vals, index_base, is_sorted);

    return matrix;
}


} //end namespace

