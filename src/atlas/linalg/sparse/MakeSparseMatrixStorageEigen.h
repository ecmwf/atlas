/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <any>
#include <memory>
#include <algorithm>
#include <type_traits>

#include "atlas/library/defines.h"
#if ATLAS_HAVE_EIGEN
ATLAS_SUPPRESS_WARNINGS_PUSH
ATLAS_SUPPRESS_WARNINGS_INTEGER_SIGN_CHANGE
ATLAS_SUPPRESS_WARNINGS_CODE_IS_UNREACHABLE
ATLAS_SUPPRESS_WARNINGS_UNUSED_BUT_SET_VARIABLE
#include <Eigen/Sparse>
ATLAS_SUPPRESS_WARNINGS_POP
#endif

#include "atlas/array.h"
#include "atlas/linalg/sparse/SparseMatrixStorage.h"

namespace atlas {
namespace linalg {

//----------------------------------------------------------------------------------------------------------------------

#if ATLAS_HAVE_EIGEN

template<typename OutputT, typename InputT>
void host_copy_eigen (const InputT* input_data, atlas::array::Array& output) {
        auto size = output.size();
        OutputT* output_data = output.host_data<OutputT>();
        std::copy( input_data, input_data + size, output_data );
}

template <typename Value, typename Index>
SparseMatrixStorage make_sparse_matrix_storage(Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>&& m) {

    std::size_t rows = m.rows();
    std::size_t cols = m.cols();
    std::size_t nnz  = m.nonZeros();
    auto* outer = atlas::array::Array::wrap(const_cast<Index*>(m.outerIndexPtr()), atlas::array::make_shape(rows+1));
    auto* inner = atlas::array::Array::wrap(const_cast<Index*>(m.innerIndexPtr()), atlas::array::make_shape(nnz));
    auto* value = atlas::array::Array::wrap(const_cast<Value*>(m.valuePtr()),      atlas::array::make_shape(nnz));

    using EigenMatrix = Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>;
    auto m_ptr = std::make_shared<EigenMatrix>();
    m_ptr->swap(m);
    auto storage = std::make_any<std::shared_ptr<EigenMatrix>>(std::move(m_ptr));

    return SparseMatrixStorage::make( rows, cols, nnz,
        std::unique_ptr<atlas::array::Array>(value),
        std::unique_ptr<atlas::array::Array>(inner),
        std::unique_ptr<atlas::array::Array>(outer),
        std::move(storage)
    );
}

template< typename value_type, typename index_type = eckit::linalg::Index, typename Value, typename Index,
          typename = std::enable_if_t < !std::is_same_v<value_type,Value>>>
SparseMatrixStorage make_sparse_matrix_storage(Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>&& m) {

    if (std::is_same_v<Value, value_type> && std::is_same_v<Index, index_type>) {
        return make_sparse_matrix_storage(std::move(m));
    }
    std::size_t rows = m.rows();
    std::size_t cols = m.cols();
    std::size_t nnz  = m.nonZeros();
    auto* outer = atlas::array::Array::create<index_type>(rows+1);
    auto* inner = atlas::array::Array::create<index_type>(nnz);
    auto* value = atlas::array::Array::create<value_type>(nnz);

    host_copy_eigen<index_type>(m.outerIndexPtr(), *outer);
    host_copy_eigen<index_type>(m.innerIndexPtr(), *inner);
    host_copy_eigen<value_type>(m.valuePtr(),      *value);

    return SparseMatrixStorage::make( rows, cols, nnz,
        std::unique_ptr<atlas::array::Array>(value),
        std::unique_ptr<atlas::array::Array>(inner),
        std::unique_ptr<atlas::array::Array>(outer),
        std::any()
    );
}


template <typename Value, typename Index>
SparseMatrixStorage make_sparse_matrix_storage(const Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>& m) {
    std::size_t rows = m.rows();
    std::size_t cols = m.cols();
    std::size_t nnz  = m.nonZeros();
    auto* outer = atlas::array::Array::create<Index>(rows+1);
    auto* inner = atlas::array::Array::create<Index>(nnz);
    auto* value = atlas::array::Array::create<Value>(nnz);

    host_copy_eigen<Index>(m.outerIndexPtr(), *outer);
    host_copy_eigen<Index>(m.innerIndexPtr(), *inner);
    host_copy_eigen<Value>(m.valuePtr(),      *value);

    return SparseMatrixStorage::make( rows, cols, nnz,
        std::unique_ptr<atlas::array::Array>(value),
        std::unique_ptr<atlas::array::Array>(inner),
        std::unique_ptr<atlas::array::Array>(outer),
        std::any()
    );
}

template< typename value_type, typename index_type = eckit::linalg::Index, typename Value, typename Index,
          typename = std::enable_if_t < !std::is_same_v<value_type,Value>>>
SparseMatrixStorage make_sparse_matrix_storage(const Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>& m) {

    std::size_t rows = m.rows();
    std::size_t cols = m.cols();
    std::size_t nnz  = m.nonZeros();
    auto* outer = atlas::array::Array::create<index_type>(rows+1);
    auto* inner = atlas::array::Array::create<index_type>(nnz);
    auto* value = atlas::array::Array::create<value_type>(nnz);
    host_copy_eigen<index_type>(m.outerIndexPtr(), *outer);
    host_copy_eigen<index_type>(m.innerIndexPtr(), *inner);
    host_copy_eigen<value_type>(m.valuePtr(),      *value);

    return SparseMatrixStorage::make( rows, cols, nnz,
        std::unique_ptr<atlas::array::Array>(value),
        std::unique_ptr<atlas::array::Array>(inner),
        std::unique_ptr<atlas::array::Array>(outer),
        std::any()
    );
}

#endif

//----------------------------------------------------------------------------------------------------------------------

}  // namespace linalg
}  // namespace atlas