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
    static_assert(std::is_integral_v<Index>, "Index must be an integral type");

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

namespace detail {
// SFINAE types/variables in place of C++20 concepts.
template <typename T>
struct is_triplet : std::false_type {};

template <typename Value, typename Index>
struct is_triplet<Triplet<Value, Index>> : std::true_type {};

template <typename Iter>
constexpr bool is_triplet_iterator = is_triplet<typename std::iterator_traits<Iter>::value_type>::value;

template <typename Iter>
constexpr bool is_random_access_iterator =
    std::is_base_of_v<typename std::iterator_traits<Iter>::iterator_category, std::random_access_iterator_tag>;

template <typename Iter>
constexpr bool is_mutable_iterator =
    !std::is_const_v<typename std::remove_reference_t<typename std::iterator_traits<Iter>::reference>>;

template <typename Iter>
constexpr bool is_sortable_iterator = is_random_access_iterator<Iter> && is_mutable_iterator<Iter>;
}  // namespace detail

/// @brief Construct a SparseMatrixStorage from a range of triplets.
template <typename Iter>
std::enable_if_t<detail::is_triplet_iterator<Iter>, SparseMatrixStorage> make_sparse_matrix_storage_from_triplets(
    std::size_t n_rows, std::size_t n_cols, Iter triplets_begin, Iter triplets_end) {
    using TripletType = typename std::iterator_traits<Iter>::value_type;
    using Index       = decltype(std::declval<TripletType>().row());
    using Value       = decltype(std::declval<TripletType>().value());

    const std::size_t n_non_zero = std::distance(triplets_begin, triplets_end);
    auto outer_array             = std::unique_ptr<array::Array>(array::Array::create<Index>(n_rows + 1));
    auto inner_array             = std::unique_ptr<array::Array>(array::Array::create<Index>(n_non_zero));
    auto values_array            = std::unique_ptr<array::Array>(array::Array::create<Value>(n_non_zero));
    auto outer_view              = array::make_view<Index, 1>(*outer_array);
    auto inner_view              = array::make_view<Index, 1>(*inner_array);
    auto values_view             = array::make_view<Value, 1>(*values_array);

    std::size_t index            = 0;
    Iter triplet_iter            = triplets_begin;
    TripletType previous_triplet = *triplet_iter;
    for (std::size_t row = 0; row < n_rows; ++row) {
        outer_view(row) = index;
        for (; triplet_iter != triplets_end && triplet_iter->row() == static_cast<Index>(row);
             ++index, ++triplet_iter) {
            ATLAS_ASSERT(!(*triplet_iter < previous_triplet), "Triplet range must be sorted.");
            ATLAS_ASSERT(static_cast<std::size_t>(triplet_iter->col()) < n_cols, "Triplet column index out of bounds.");
            inner_view(index)  = triplet_iter->col();
            values_view(index) = triplet_iter->value();
            previous_triplet   = *triplet_iter;
        }
    }
    ATLAS_ASSERT(index == n_non_zero, "Triplet row index out of bounds.");
    outer_view(n_rows) = n_non_zero;

    return SparseMatrixStorage::make(n_rows, n_cols, n_non_zero, std::move(values_array), std::move(inner_array),
                                     std::move(outer_array), std::any());
}

/// @brief Construct a SparseMatrixStorage from a random-access range of triplets, sorting if necessary.
template <typename Iter>
std::enable_if_t<detail::is_triplet_iterator<Iter> && detail::is_sortable_iterator<Iter>, SparseMatrixStorage>
make_sparse_matrix_storage_from_triplets(std::size_t n_rows, std::size_t n_cols, std::size_t n_non_zero,
                                         Iter triplets_begin, bool is_sorted = false) {
    if (!is_sorted) {
        std::sort(triplets_begin, triplets_begin + n_non_zero);
    }
    return make_sparse_matrix_storage_from_triplets(n_rows, n_cols, triplets_begin, triplets_begin + n_non_zero);
}

/// @brief Construct a SparseMatrixStorage from a vector of triplets, sorting if necessary.
template <typename Value, typename Index>
SparseMatrixStorage make_sparse_matrix_storage_from_triplets(std::size_t n_rows, std::size_t n_cols,
                                                             std::vector<Triplet<Value, Index>>& triplets,
                                                             bool is_sorted = false) {
    return make_sparse_matrix_storage_from_triplets(n_rows, n_cols, triplets.size(), triplets.begin(), is_sorted);
}

/// @brief For-each iteration over all non-zero triplets in row.
template <typename Value, typename Index, typename Functor>
std::enable_if_t<std::is_invocable_v<Functor, Index, Index, Value>> sparse_matrix_for_each_row(
    std::size_t row, const SparseMatrixView<Value, Index>& matrix, Functor&& functor) {
    const Index* outer  = matrix.outer();
    const Index* inner  = matrix.inner();
    const Value* values = matrix.value();

    for (Index index = outer[row]; index < outer[row + 1]; ++index) {
        const Index column = inner[index];
        const Value value  = values[index];
        functor(static_cast<Index>(row), column, value);
    }
}

/// @brief For-each iteration over all non-zero triplets in matrix.
template <typename Value, typename Index, typename Functor>
std::enable_if_t<std::is_invocable_v<Functor, Index, Index, Value>> sparse_matrix_for_each(
    const SparseMatrixView<Value, Index>& matrix, Functor&& functor) {
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        sparse_matrix_for_each_row(row, matrix, functor);
    }
}

}  // namespace atlas::linalg
