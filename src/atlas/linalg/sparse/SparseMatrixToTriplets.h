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
#pragma once

#include <vector>
#include <numeric>

#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/linalg/sparse/SparseMatrixView.h"
#include "atlas/util/DataType.h"
#include "atlas/runtime/Trace.h"

namespace atlas::linalg {

template<typename ViewIndex, typename ViewValue, typename Index, typename Value>
void sparse_matrix_to_rows_columns_values(const linalg::SparseMatrixView<ViewValue, ViewIndex>& mat, Index rows[], Index columns[], Value values[]) {
    const auto outer  = mat.outer();
    const auto index  = mat.inner();
    const auto vals   = mat.value();

    std::size_t idx = 0;
    for(std::size_t r = 0; r < mat.rows(); ++r) {
        for (auto c = outer[r]; c < outer[r + 1]; ++c) {
            auto col = index[c];
            rows[idx] = r;
            columns[idx] = col;
            values[idx] = vals[c];
            ++idx;
        }
    }
}


template<typename ViewIndex, typename ViewValue, typename Index, typename Value>
void sparse_matrix_to_rows_columns_values(const linalg::SparseMatrixView<ViewValue, ViewIndex>& mat, std::vector<Index>& rows, std::vector<Index>& columns, std::vector<Value>& values) {
    rows.resize(mat.nnz());
    columns.resize(mat.nnz());
    values.resize(mat.nnz());
    sparse_matrix_to_rows_columns_values(mat, rows.data(), columns.data(), values.data());
}


template<typename Index, typename Value>
void sparse_matrix_to_rows_columns_values(const SparseMatrixStorage& mat, std::vector<Index>& rows, std::vector<Index>& columns, std::vector<Value>& values) {
    if (mat.inner().datatype().kind() != DataType::kind<Index>() || mat.outer().datatype().kind() != DataType::kind<Index>() ) {
        ATLAS_NOTIMPLEMENTED;
    }
    switch (mat.value().datatype().kind() ) {
        case DataType::kind<double>() : {
            sparse_matrix_to_rows_columns_values( atlas::linalg::make_host_view<double, Index>(mat), rows, columns, values );
            return;
        }
        case DataType::kind<float>() : {
            sparse_matrix_to_rows_columns_values( atlas::linalg::make_host_view<float, Index>(mat), rows, columns, values );
            return;
        }
        default:
            ATLAS_NOTIMPLEMENTED;
    }
}

template<typename SparseMatrixValue = eckit::linalg::Scalar, typename SparseMatrixIndex = eckit::linalg::Index, typename Value, typename Index, typename IndexBase>
SparseMatrixStorage make_sparse_matrix_storage_from_rows_columns_values(std::size_t nr, std::size_t nc, std::size_t nnz, const Index rows[], const Index cols[], const Value vals[], const IndexBase index_base = 0, bool is_sorted = true) {
    std::unique_ptr<array::Array> array_value(array::Array::create<SparseMatrixValue>(nnz));
    std::unique_ptr<array::Array> array_inner(array::Array::create<SparseMatrixIndex>(nnz));
    std::unique_ptr<array::Array> array_outer(array::Array::create<SparseMatrixIndex>(nr+1));
    auto* value = array_value->host_data<SparseMatrixValue>();
    auto* inner = array_inner->host_data<SparseMatrixIndex>();
    auto* outer = array_outer->host_data<SparseMatrixIndex>();
    std::fill(outer, outer + nr + 1, 0);

    Index base = index_base;

    if (not is_sorted) {
        std::vector<size_t> sorted_index(nnz);
        std::iota(sorted_index.begin(), sorted_index.end(), 0);
        std::sort(sorted_index.begin(), sorted_index.end(), [&](auto i, auto j) {
            return rows[i] != rows[j] ? rows[i] < rows[j] :
                   cols[i] != cols[j] ? cols[i] < cols[j] :
                   i < j;
        });

        for (size_t n = 0; n < nnz; ++n) {
            auto r = rows[sorted_index[n]] - base;
            auto c = cols[sorted_index[n]] - base;
            outer[r + 1]++;
            inner[n] = c;
        }
        for (size_t n = 0; n < nnz; ++n) {
            value[n] = vals[sorted_index[n]];
        }
    }
    else {
        for (size_t n = 0; n < nnz; ++n) {
            auto r = rows[n] - base;
            auto c = cols[n] - base;
            outer[r + 1]++;
            inner[n] = c;
        }
        for (size_t n = 0; n < nnz; ++n) {
            value[n] = vals[n];
        }

    }
    for (size_t r = 0; r < nr; ++r) {
        outer[r + 1] += outer[r];
    }
    ATLAS_ASSERT(outer[0] == 0);
    ATLAS_ASSERT(outer[nr] == nnz);
    return SparseMatrixStorage::make(nr, nc, nnz, std::move(array_value), std::move(array_inner), std::move(array_outer), std::any());
}

template<typename SparseMatrixValue = eckit::linalg::Scalar, typename SparseMatrixIndex = eckit::linalg::Index, typename Value, typename Index, typename IndexBase>
SparseMatrixStorage make_sparse_matrix_storage_from_rows_columns_values(std::size_t nr, std::size_t nc, const std::vector<Index>& rows, const std::vector<Index>& cols, const std::vector<Value>& vals, const IndexBase index_base = 0, bool is_sorted = true) {
    ATLAS_TRACE("make_sparse_matrix_storage_from_rows_columns_values partition");
    std::size_t nnz = vals.size();
    ATLAS_ASSERT(rows.size() == nnz);
    ATLAS_ASSERT(cols.size() == nnz);
    return make_sparse_matrix_storage_from_rows_columns_values<SparseMatrixValue, SparseMatrixIndex>(nr, nc, nnz, rows.data(), cols.data(), vals.data(), index_base, is_sorted);
}


} // namespace atlas::linalg
