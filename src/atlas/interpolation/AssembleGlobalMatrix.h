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

} // namespace atlas::interpolation
