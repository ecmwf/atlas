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
#include <vector>

#include "eckit/linalg/types.h"

#include "atlas/linalg/sparse/SparseMatrixStorage.h"

namespace atlas::grid { class Distribution; }

namespace atlas::interpolation {

    atlas::linalg::SparseMatrixStorage assemble_global_matrix(size_t rows, size_t cols, const Interpolation& interpolation, int mpi_root = 0);
    atlas::linalg::SparseMatrixStorage assemble_global_matrix(const Interpolation& interpolation, int mpi_root = 0);

    atlas::linalg::SparseMatrixStorage distribute_global_matrix(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs, const linalg::SparseMatrixStorage&, int mpi_root = 0);

    /// Send each row of the global matrix (only needed on mpi_root) as (row, col, value) triplets
    /// with 0-based global indices to the task tgt_partition[row]. tgt_partition (size = rows)
    /// is only accessed on mpi_root.
    void distribute_global_matrix_as_triplets(const int* tgt_partition, const linalg::SparseMatrixStorage& gmatrix,
        std::vector<eckit::linalg::Index>& rows, std::vector<eckit::linalg::Index>& cols,
        std::vector<eckit::linalg::Scalar>& vals, int mpi_root = 0);

    void distribute_global_matrix_as_triplets(const grid::Distribution& tgt_distribution, const linalg::SparseMatrixStorage& gmatrix,
        std::vector<eckit::linalg::Index>& rows, std::vector<eckit::linalg::Index>& cols,
        std::vector<eckit::linalg::Scalar>& vals, int mpi_root = 0);

    atlas::linalg::SparseMatrixStorage distribute_global_matrix(const grid::Distribution& tgt_distribution, const FunctionSpace& src_fs, const FunctionSpace& tgt_fs, const linalg::SparseMatrixStorage&, int mpi_root = 0);

} // namespace atlas::interpolation
