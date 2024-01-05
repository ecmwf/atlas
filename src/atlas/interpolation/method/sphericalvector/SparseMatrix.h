/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

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
  using Triplets = std::vector<Eigen::Triplet<Value>>;

  SparseMatrix(Index nRows, Index nCols, const Triplets& triplets)
      : eigenMatrix_(nRows, nCols) {
    eigenMatrix_.setFromTriplets(triplets.begin(), triplets.end());
  }

  Size nonZeros() const { return eigenMatrix_.nonZeros(); }
  Size rows() const { return eigenMatrix_.rows(); }
  Size cols() const { return eigenMatrix_.cols(); }
  const Index* outer() { return eigenMatrix_.outerIndexPtr(); }
  const Index* inner() { return eigenMatrix_.innerIndexPtr(); }
  const Value* data() { return eigenMatrix_.valuePtr(); }
  SparseMatrix<Value> adjoint() {
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
  class Triplet {
   public:
    template <typename... Args>
    Triplet(const Args&... args) {}
  };
  using Index = int;
  using Size = long int;
  using Triplets = std::vector<Triplet>;

  template <typename... Args>
  SparseMatrix(const Args&... args) {
    ATLAS_THROW_EXCEPTION("Atlas has been compiled without Eigen");
  }
  constexpr Size nonZeros() const { return 0; }
  constexpr Size rows() const { return 0; }
  constexpr Size cols() const { return 0; }
  constexpr const Index* outer() { return nullptr; }
  constexpr const Index* inner() { return nullptr; }
  constexpr const Value* data() { return nullptr; }
  SparseMatrix<Value> adjoint() {
    return SparseMatrix<Value>(0, 0, Triplets{});
  }
};
#endif

} // namespace detail
} // namespace method
} // namespace interpolation
} // namespace atlas
