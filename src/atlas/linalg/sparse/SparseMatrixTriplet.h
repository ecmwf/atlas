/*
 * (C) Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <cstring>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/linalg/sparse/SparseMatrixView.h"


namespace atlas::linalg {

/// @brief A triplet which represents a non-zero entry in a sparse matrix.
template <typename Value, typename Index>
class Triplet {
public:
    Triplet() = default;
    Triplet(Index row, Index column, Value value): row_{row}, column_{column}, value_{value} {}
    Index row() const { return row_; }
    Index col() const { return column_; }
    Value value() const { return value_; }

    /// @brief Less-than operator for sorting triplets.
    bool operator<(const Triplet& other) const { return std::tie(row_, column_) < std::tie(other.row_, other.column_); }

private:
    Index row_{};
    Index column_{};
    Value value_{};
};

template <typename Value, typename Index>
SparseMatrixStorage make_sparse_matrix_storage_from_triplets(Index n_rows, Index n_cols,
                                                             std::vector<Triplet<Value, Index>>&& triplets) {
    const auto n_non_zero = triplets.size();
    auto outer_array      = std::unique_ptr<array::Array>(array::Array::create<Index>(n_rows + 1));
    auto inner_array      = std::unique_ptr<array::Array>(array::Array::create<Index>(triplets.size()));
    auto values_array     = std::unique_ptr<array::Array>(array::Array::create<Value>(triplets.size()));
    auto outer_view       = array::make_view<Index, 1>(*outer_array);
    auto inner_view       = array::make_view<Index, 1>(*inner_array);
    auto values_view      = array::make_view<Value, 1>(*values_array);

    std::sort(triplets.begin(), triplets.end());
    ATLAS_ASSERT(triplets.back().row() < n_rows);

    std::size_t index = 0;
    for (Index row = 0; row < n_rows; ++row) {
        outer_view(row) = index;
        while (index < triplets.size() && triplets[index].row() == row) {
            ATLAS_ASSERT(triplets[index].col() < n_cols);
            inner_view(index)  = triplets[index].col();
            values_view(index) = triplets[index].value();
            ++index;
        }
    }
    ATLAS_ASSERT(index == n_non_zero);

    outer_view(n_rows) = n_non_zero;
    triplets.clear(); // Leave vector in empty state after move.

    return SparseMatrixStorage::make(n_rows, n_cols, n_non_zero, std::move(values_array), std::move(inner_array),
                                     std::move(outer_array), std::any());
}

// For-each iteration over all non-zero elements in row.
template <typename Value, typename Index, typename Functor>
std::enable_if_t<std::is_invocable_v<Functor, const Triplet<Value, Index>>> sparse_matrix_for_each_triplet(
    std::size_t row, const SparseMatrixView<Value, Index>& matrix, Functor&& functor) {
    const Index* outer  = matrix.outer();
    const Index* inner  = matrix.inner();
    const Value* values = matrix.value();

    for (auto index = outer[row]; index < outer[row + 1]; ++index) {
        Index column = inner[index];
        Value value  = values[index];
        functor(Triplet{static_cast<Index>(row), column, value});
    }
}

// For-each iteration over all non-zero elements in matrix.
template <typename Value, typename Index, typename Functor>
std::enable_if_t<std::is_invocable_v<Functor, const Triplet<Value, Index>>> sparse_matrix_for_each_triplet(
    const SparseMatrixView<Value, Index>& matrix, Functor&& functor) {
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        sparse_matrix_for_each_triplet(row, matrix, functor);
    }
}


}  // namespace atlas::linalg