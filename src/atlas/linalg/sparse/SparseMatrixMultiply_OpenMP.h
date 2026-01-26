/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/linalg/sparse/SparseMatrixMultiply.h"

namespace atlas {
namespace linalg {
namespace sparse {


template <typename MatrixValue, typename SourceValue, typename TargetValue>
struct SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 1, MatrixValue, SourceValue, TargetValue> {
    static void multiply(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 1>& src, View<TargetValue, 1>& tgt,
                      const Configuration&);
    static void multiply_add(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 1>& src, View<TargetValue, 1>& tgt,
                      const Configuration&);
};

template <typename MatrixValue, typename SourceValue, typename TargetValue>
struct SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 2, MatrixValue, SourceValue, TargetValue> {
    static void multiply(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 2>& src, View<TargetValue, 2>& tgt,
                      const Configuration&);
    static void multiply_add(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 2>& src, View<TargetValue, 2>& tgt,
                      const Configuration&);
};

template <typename MatrixValue, typename SourceValue, typename TargetValue>
struct SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 3, MatrixValue, SourceValue, TargetValue> {
    static void multiply(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 3>& src, View<TargetValue, 3>& tgt,
                      const Configuration&);
    static void multiply_add(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 3>& src, View<TargetValue, 3>& tgt,
                      const Configuration&);
};

template <typename MatrixValue, typename SourceValue, typename TargetValue>
struct SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 1, MatrixValue, SourceValue, TargetValue> {
    static void multiply(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 1>& src, View<TargetValue, 1>& tgt,
                      const Configuration&);
    static void multiply_add(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 1>& src, View<TargetValue, 1>& tgt,
                      const Configuration&);
};

template <typename MatrixValue, typename SourceValue, typename TargetValue>
struct SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 2, MatrixValue, SourceValue, TargetValue> {
    static void multiply(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 2>& src, View<TargetValue, 2>& tgt,
                      const Configuration&);
    static void multiply_add(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 2>& src, View<TargetValue, 2>& tgt,
                      const Configuration&);
};

template <typename MatrixValue, typename SourceValue, typename TargetValue>
struct SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 3, MatrixValue, SourceValue, TargetValue> {
    static void multiply(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 3>& src, View<TargetValue, 3>& tgt,
                      const Configuration&);
    static void multiply_add(const SparseMatrixView<MatrixValue>& W, const View<SourceValue, 3>& src, View<TargetValue, 3>& tgt,
                      const Configuration&);
};

}  // namespace sparse
}  // namespace linalg
}  // namespace atlas
