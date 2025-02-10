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

#include "eckit/linalg/SparseMatrix.h"

#include "atlas/array.h"
#include "atlas/linalg/sparse/SparseMatrixStorage.h"

namespace atlas {
namespace linalg {

//----------------------------------------------------------------------------------------------------------------------

template<typename OutputT, typename InputT>
void host_copy_eckit (const InputT* input_data, array::Array& output) {
        auto size = output.size();
        OutputT* output_data = output.host_data<OutputT>();
        std::copy( input_data, input_data + size, output_data );
}

inline SparseMatrixStorage make_sparse_matrix_storage(eckit::linalg::SparseMatrix&& m) {
    using Index = eckit::linalg::Index;
    using Value = eckit::linalg::Scalar;

    std::size_t rows = m.rows();
    std::size_t cols = m.cols();
    std::size_t nnz  = m.nonZeros();
    auto* outer = array::Array::wrap(const_cast<Index*>(m.outer()), array::make_shape(rows+1));
    auto* inner = array::Array::wrap(const_cast<Index*>(m.inner()), array::make_shape(nnz));
    auto* value = array::Array::wrap(const_cast<Value*>(m.data()),  array::make_shape(nnz));

    // We now move the eckit::linalg::SparseMatrix into a generic storage so
    //   the wrapped array data does not go out of scope
    auto storage = std::make_any<eckit::linalg::SparseMatrix>(std::move(m));

    return SparseMatrixStorage::make( rows, cols, nnz,
        std::unique_ptr<array::Array>(value),
        std::unique_ptr<array::Array>(inner),
        std::unique_ptr<array::Array>(outer),
        std::move(storage)
    );
}

template< typename value_type, typename index_type = eckit::linalg::Index>
SparseMatrixStorage make_sparse_matrix_storage(eckit::linalg::SparseMatrix&& m) {
    if (std::is_same_v<eckit::linalg::Scalar, value_type> && std::is_same_v<eckit::linalg::Index, index_type>) {
        return make_sparse_matrix_storage(std::move(m));
    }
    std::size_t rows = m.rows();
    std::size_t cols = m.cols();
    std::size_t nnz  = m.nonZeros();
    auto* outer = array::Array::create<index_type>(rows+1);
    auto* inner = array::Array::create<index_type>(nnz);
    auto* value = array::Array::create<value_type>(nnz);

    host_copy_eckit<index_type>(m.outer(), *outer);
    host_copy_eckit<index_type>(m.inner(), *inner);
    host_copy_eckit<value_type>(m.data(),  *value);

    return SparseMatrixStorage::make( rows, cols, nnz,
        std::unique_ptr<array::Array>(value),
        std::unique_ptr<array::Array>(inner),
        std::unique_ptr<array::Array>(outer),
        std::any()
    );
}


inline SparseMatrixStorage make_sparse_matrix_storage(const eckit::linalg::SparseMatrix& m) {
    using Index = eckit::linalg::Index;
    using Value = eckit::linalg::Scalar;

    std::size_t rows = m.rows();
    std::size_t cols = m.cols();
    std::size_t nnz  = m.nonZeros();
    auto* outer = array::Array::create<Index>(rows+1);
    auto* inner = array::Array::create<Index>(nnz);
    auto* value = array::Array::create<Value>(nnz);

    host_copy_eckit<Index>(m.outer(), *outer);
    host_copy_eckit<Index>(m.inner(), *inner);
    host_copy_eckit<Value>(m.data(),  *value);

    return SparseMatrixStorage::make( rows, cols, nnz,
        std::unique_ptr<array::Array>(value),
        std::unique_ptr<array::Array>(inner),
        std::unique_ptr<array::Array>(outer),
        std::any()
    );
}

template< typename value_type, typename index_type = eckit::linalg::Index>
SparseMatrixStorage make_sparse_matrix_storage(const eckit::linalg::SparseMatrix& m) {

    std::size_t rows = m.rows();
    std::size_t cols = m.cols();
    std::size_t nnz  = m.nonZeros();
    auto* outer = array::Array::create<index_type>(rows+1);
    auto* inner = array::Array::create<index_type>(nnz);
    auto* value = array::Array::create<value_type>(nnz);
    host_copy_eckit<index_type>(m.outer(), *outer);
    host_copy_eckit<index_type>(m.inner(), *inner);
    host_copy_eckit<value_type>(m.data(),  *value);

    return SparseMatrixStorage::make( rows, cols, nnz,
        std::unique_ptr<array::Array>(value),
        std::unique_ptr<array::Array>(inner),
        std::unique_ptr<array::Array>(outer),
        std::any()
    );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace linalg
}  // namespace atlas