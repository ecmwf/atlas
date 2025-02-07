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

namespace atlas::interpolation {

    atlas::linalg::SparseMatrixStorage assemble_global_matrix(const Interpolation& interpolation, int mpi_root = 0);

    atlas::linalg::SparseMatrixStorage distribute_global_matrix(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs, const linalg::SparseMatrixStorage&, int mpi_root = 0);

} // namespace atlas::interpolation
