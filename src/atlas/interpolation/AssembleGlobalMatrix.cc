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

#include "eckit/linalg/types.h"

#include "atlas/array.h"
#include "atlas/linalg/sparse/SparseMatrixToTriplets.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/interpolation/Cache.h"

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

            const auto src_remote_index = array::make_indexview<idx_t, 1>(src_fs.remote_index());
            idx_t one_if_structuredcolumns = (src_fs.type() == "StructuredColumns" ? 1 : 0);
            auto src_ridx = [&](auto idx) -> idx_t {
                return src_remote_index(idx) + one_if_structuredcolumns;
            };

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
        mpi_comm.reduceInPlace(max_gidx, eckit::mpi::max(), mpi_root);
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


linalg::SparseMatrixStorage distribute_global_matrix(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs, const linalg::SparseMatrixStorage& gmatrix, int mpi_root) {
    const auto src_global_index = array::make_view<gidx_t, 1>(src_fs.global_index());
    const auto tgt_global_index = array::make_view<gidx_t, 1>(tgt_fs.global_index());
    const auto src_part         = array::make_view<int, 1>(src_fs.partition());
    const auto tgt_part         = array::make_view<int, 1>(tgt_fs.partition());

    Field field_tgt_part_glb = tgt_fs.createField(tgt_fs.partition(), option::global(mpi_root) );
    ATLAS_DEBUG_VAR(field_tgt_part_glb.size());
    ATLAS_DEBUG_VAR(tgt_fs.partition().size());
    tgt_fs.gather(tgt_fs.partition(), field_tgt_part_glb);
    const auto tgt_part_glb = array::make_view<int,1>(field_tgt_part_glb);
    
    auto& mpi_comm = mpi::comm();
    auto mpi_size  = mpi_comm.size();
    auto mpi_rank  = mpi_comm.rank();

    std::cout << mpi_rank << " src nb_part " << src_fs.nb_parts() << ", size " << src_fs.partition().size() << std::endl;
    std::cout << mpi_rank << " tgt nb_part " << tgt_fs.nb_parts() << ", size " << tgt_fs.partition().size() << std::endl;

    const auto src_remote_index = array::make_indexview<idx_t, 1>(src_fs.remote_index());
    idx_t one_if_structuredcolumns = (src_fs.type() == "StructuredColumns" ? 1 : 0);
    auto src_ridx = [&](auto idx) -> idx_t {
        return src_remote_index(idx) + one_if_structuredcolumns;
    };

    using Scalar = eckit::linalg::Scalar;
    using Index = eckit::linalg::Index;
    std::vector<Index> rows, cols;
    std::vector<Scalar> vals;

    // compute how many nnz-entries each task gets
    size_t nnz_loc = 0;
    int mpi_tag = 0;
    if (mpi_rank == mpi_root) {
        auto global_matrix = linalg::make_host_view<Scalar, Index>(gmatrix);
        
        const auto outer  = global_matrix.outer();
        const auto inner  = global_matrix.inner();
        const auto value  = global_matrix.value();

        std::vector<std::size_t> nnz_per_task(mpi_size);
        for(std::size_t r = 0; r < global_matrix.rows(); ++r) {
            //std::cout << "r, tgt_part: " << r << ", part " << tgt_part_glb(r) << ", cols " << outer[r+1] - outer[r] << std::endl;
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

    linalg::SparseMatrixStorage matrix;

    mpi_comm.barrier();

    return matrix;
/*
    size_t idx = 0;
    for(std::size_t r = 0; r < global_matrix.rows(); ++r) {
        for (auto c = outer[r]; c < outer[r + 1]; ++c) {
            auto col = inner[c];
            rows[idx] = r;
            columns[idx] = col;
            values[idx] = vals[c];
            ++idx;
        }
    }

    for( row in global_matrix)
    
    //linalg::SparseMatrixStorage local_matrix(src_fs.size(), tgt_fs.size());
    // row_gl, col_gl -> (row, spart), (col, tpart)
    //
    for (int i = 0; i < tgt_fs.size(); ++i) {
    }

    Config config;
    config.setString("type", "finite-element");

    interpolation::MatrixCache cache(std::move(local_matrix));
    // now don't use local_matrix anymore.

    auto interp = Interpolation(config, src_fs, tgt_fs, cache);
*/
}


} //end namespace

