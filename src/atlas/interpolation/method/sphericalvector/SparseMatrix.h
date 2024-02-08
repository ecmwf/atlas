/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "atlas/library/defines.h"
#if ATLAS_HAVE_EIGEN
#include <Eigen/Sparse>
#endif

#include "atlas/runtime/Exception.h"

namespace atlas {
namespace interpolation {
namespace method {
namespace detail {

#if ATLAS_HAVE_EIGEN
/// @brief   Wrapper class for Eigen sparse matrix
///
/// @details Adapts the Eigen sparse matrix interface to be more in line with
///          eckit::linalg::SparseMatrix. Also allows preprocessor disabling of
///          class is Eigen library is not present.
template <typename Value>
class SparseMatrix {
  using EigenMatrix = Eigen::SparseMatrix<Value, Eigen::RowMajor>;

 public:
  using Index = typename EigenMatrix::StorageIndex;
  using Size = typename EigenMatrix::Index;
  using Triplet = Eigen::Triplet<Value>;
  using Triplets = std::vector<Triplet>;
  using RowIter = typename EigenMatrix::InnerIterator;

  SparseMatrix(Index nRows, Index nCols, const Triplets& triplets)
      : eigenMatrix_(nRows, nCols) {
    eigenMatrix_.setFromTriplets(triplets.begin(), triplets.end());
  }

  Size nonZeros() const { return eigenMatrix_.nonZeros(); }
  Size rows() const { return eigenMatrix_.rows(); }
  Size cols() const { return eigenMatrix_.cols(); }
  RowIter rowIter(Size rowIndex) const {
    return RowIter(eigenMatrix_, rowIndex);
  }
  SparseMatrix<Value> adjoint() const {
    return SparseMatrix(eigenMatrix_.adjoint().eval());
  }

 private:
  SparseMatrix(EigenMatrix&& eigenMatrixAdjoint)
      : eigenMatrix_{std::move(eigenMatrixAdjoint)} {}
  EigenMatrix eigenMatrix_{};
};
#else

template <typename Value>
class SparseMatrix {
 public:
  using Index = int;
  using Size = long int;

  class Triplet {
   public:
    template <typename... Args>
    constexpr Triplet(const Args&... args) {}
  };
  using Triplets = std::vector<Triplet>;

  class RowIter {
   public:
    template <typename... Args>
    constexpr RowIter(const Args&... args) {}
    constexpr Index row() const { return Index{}; }
    constexpr Index col() const { return Index{}; }
    constexpr Value value() const { return Value{}; }
    constexpr operator bool() const { return false; }
    constexpr RowIter& operator++() { return *this; }
  };

  template <typename... Args>
  SparseMatrix(const Args&... args) {
    throw_Exception("Atlas has been compiled without Eigen", Here());
  }
  constexpr Size nonZeros() const { return Size{}; }
  constexpr Size rows() const { return Size{}; }
  constexpr Size cols() const { return Size{}; }
  constexpr RowIter rowIter(Size rowIndex) const { return RowIter{}; }
  constexpr SparseMatrix<Value> adjoint() const {
    return SparseMatrix<Value>{};
  }
};
#endif

}  // namespace detail
}  // namespace method
}  // namespace interpolation
}  // namespace atlas
