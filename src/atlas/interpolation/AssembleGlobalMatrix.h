/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/interpolation/Interpolation.h"
#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/grid/Distribution.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"

namespace atlas::grid { class Distribution; }

namespace atlas::interpolation {

    linalg::SparseMatrixStorage assemble_global_matrix(const Interpolation& interpolation, int mpi_root = 0);

    template <typename partition_t, typename ViewValue, typename ViewIndex, typename Value, typename Index>
    void distribute_global_matrix_as_triplets(partition_t tgt_partition, const linalg::SparseMatrixView<ViewValue,ViewIndex>& global_matrix,
    std::vector<Index>& rows, std::vector<Index>& cols, std::vector<Value>& vals, int mpi_root, std::string tgt_comm_name = "world");

    template <typename ViewValue, typename ViewIndex, typename Value, typename Index>
    void distribute_global_matrix_as_triplets( const array::Array& tgt_partition,
        const linalg::SparseMatrixView<ViewValue,ViewIndex>& global_matrix,
        std::vector<Index>& rows, std::vector<Index>& cols, std::vector<Value>& vals, int mpi_root, std::string tgt_mpi_comm = "world");

    template <typename ViewValue, typename ViewIndex, typename Value, typename Index>
    void distribute_global_matrix_as_triplets(
        const grid::Distribution& tgt_distribution, const linalg::SparseMatrixView<ViewValue,ViewIndex>& global_matrix,
        std::vector<Index>& rows, std::vector<Index>& cols, std::vector<Value>& vals, int mpi_root, std::string tgt_mpi_comm = "world");

    template<typename partition_t>
    linalg::SparseMatrixStorage distribute_global_matrix(partition_t tgt_partition, const FunctionSpace& src_fs,
        const FunctionSpace& tgt_fs, const linalg::SparseMatrixStorage& gmatrix, int mpi_root);

    linalg::SparseMatrixStorage distribute_global_matrix(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs, const linalg::SparseMatrixStorage& gmatrix, int mpi_root = 0);

    linalg::SparseMatrixStorage distribute_global_matrix(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs, const linalg::SparseMatrixStorage&, int mpi_root);

    linalg::SparseMatrixStorage distribute_global_matrix(const grid::Distribution& tgt_distribution, const FunctionSpace& src_fs, const FunctionSpace& tgt_fs, const linalg::SparseMatrixStorage&, int mpi_root = 0);

    // declarations

    template <typename partition_t, typename ViewValue, typename ViewIndex, typename Value, typename Index>
    void distribute_global_matrix_as_triplets(partition_t tgt_partition, const linalg::SparseMatrixView<ViewValue,ViewIndex>& global_matrix,
        std::vector<Index>& rows, std::vector<Index>& cols, std::vector<Value>& vals, int mpi_root, std::string tgt_comm_name) {
        ATLAS_TRACE("distribute_global_matrix_as_triplets");

        auto& mpi_comm = mpi::comm(tgt_comm_name);
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

            for(std::size_t r = 0; r < global_matrix.rows(); ++r) {
                nnz_per_task[tgt_partition[r]] += outer[r+1] - outer[r];
            }
            for (int jproc = 0; jproc < mpi_comm.size(); ++jproc) {
                if (jproc != mpi_root) {
                    mpi_comm.send(nnz_per_task.data() + jproc, 1, jproc, mpi_tag);
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
                int jproc = tgt_partition[r];
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


    template <typename ViewValue, typename ViewIndex, typename Value, typename Index>
    void distribute_global_matrix_as_triplets(
        const array::Array& tgt_partition,
        const linalg::SparseMatrixView<ViewValue,ViewIndex>& global_matrix,
        std::vector<Index>& rows, std::vector<Index>& cols, std::vector<Value>& vals, int mpi_root, std::string tgt_mpi_comm) {
        distribute_global_matrix_as_triplets(
            array::make_view<int,1>(tgt_partition).data(), global_matrix, rows, cols, vals, mpi_root, tgt_mpi_comm);
    }


    template <typename ViewValue, typename ViewIndex, typename Value, typename Index>
    void distribute_global_matrix_as_triplets(
        const grid::Distribution& tgt_distribution, const linalg::SparseMatrixView<ViewValue,ViewIndex>& global_matrix,
        std::vector<Index>& rows, std::vector<Index>& cols, std::vector<Value>& vals, int mpi_root, std::string tgt_mpi_comm) {
        struct partition_t {
            partition_t(const grid::Distribution& d) : d_(d) {}
            int operator[](idx_t i) const { return d_.partition(i); }
            const grid::Distribution& d_;
        } tgt_partition(tgt_distribution);
        distribute_global_matrix_as_triplets(tgt_partition, global_matrix, rows, cols, vals, mpi_root, tgt_mpi_comm);
    }

} // namespace atlas::interpolation
