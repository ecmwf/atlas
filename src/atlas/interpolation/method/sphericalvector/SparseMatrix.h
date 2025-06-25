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
ATLAS_SUPPRESS_WARNINGS_PUSH
ATLAS_SUPPRESS_WARNINGS_INTEGER_SIGN_CHANGE
ATLAS_SUPPRESS_WARNINGS_CODE_IS_UNREACHABLE
ATLAS_SUPPRESS_WARNINGS_UNUSED_BUT_SET_VARIABLE
#include <Eigen/Sparse>
ATLAS_SUPPRESS_WARNINGS_POP
#endif

#include "atlas/runtime/Exception.h"
#include "eckit/log/CodeLocation.h"

namespace atlas {
namespace interpolation {
namespace method {
namespace detail {

#if ATLAS_HAVE_EIGEN
/// @brief   Wrapper class for Eigen sparse matrix
///
/// @details Eigen sparse matrix wrapper. Allows preprocessor disabling of class
///          if Eigen library is not present.
template <typename Value>
class SparseMatrix {
 public:
  using Index = int;
  using EigenMatrix = Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>;
  using Triplet = Eigen::Triplet<Value, Index>;
  using Triplets = std::vector<Triplet>;
  using RowIter = typename EigenMatrix::InnerIterator;

  SparseMatrix() = default;

  SparseMatrix(Index nRows, Index nCols, const Triplets& triplets)
      : eigenMatrix_(nRows, nCols) {
    eigenMatrix_.setFromTriplets(triplets.begin(), triplets.end());
  }

  Index nonZeros() const { return eigenMatrix_.nonZeros(); }
  Index rows() const { return eigenMatrix_.rows(); }
  Index cols() const { return eigenMatrix_.cols(); }
  RowIter rowIter(Index rowIndex) const {
    return RowIter(eigenMatrix_, rowIndex);
  }
  SparseMatrix<Value> adjoint() const {
    auto adjointMatrix = SparseMatrix{};
    adjointMatrix.eigenMatrix_ = eigenMatrix_.adjoint();
    return adjointMatrix;
  }

 private:
  EigenMatrix eigenMatrix_{};
};
#else

template <typename Value>
class SparseMatrix {
 public:
  using Index = int;

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

  constexpr SparseMatrix() = default;
  template <typename... Args>
  SparseMatrix(const Args&... args) {
    throw_Exception("Atlas has been compiled without Eigen", Here());
  }
  constexpr Index nonZeros() const { return Index{}; }
  constexpr Index rows() const { return Index{}; }
  constexpr Index cols() const { return Index{}; }
  constexpr RowIter rowIter(Index rowIndex) const { return RowIter{}; }
  constexpr SparseMatrix<Value> adjoint() const {
    return SparseMatrix<Value>{};
  }
};
#endif

}  // namespace detail
}  // namespace method
}  // namespace interpolation
}  // namespace atlas
