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

#include "eckit/config/Configuration.h"

#include "atlas/linalg/Indexing.h"
#include "atlas/linalg/View.h"
#include "atlas/linalg/sparse/Backend.h"
#include "atlas/linalg/sparse/SparseMatrixView.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace linalg {

using Configuration = eckit::Configuration;

template <typename Matrix, typename SourceView, typename TargetView>
void sparse_matrix_multiply(const Matrix& matrix, const SourceView& src, TargetView& tgt);

template <typename Matrix, typename SourceView, typename TargetView>
void sparse_matrix_multiply(const Matrix& matrix, const SourceView& src, TargetView& tgt, const Configuration& config);

template <typename Matrix, typename SourceView, typename TargetView>
void sparse_matrix_multiply(const Matrix& matrix, const SourceView& src, TargetView& tgt, Indexing);

template <typename Matrix, typename SourceView, typename TargetView>
void sparse_matrix_multiply(const Matrix& matrix, const SourceView& src, TargetView& tgt, Indexing,
                            const Configuration& config);

template <typename Matrix, typename SourceView, typename TargetView>
void sparse_matrix_multiply_add(const Matrix& matrix, const SourceView& src, TargetView& tgt);

template <typename Matrix, typename SourceView, typename TargetView>
void sparse_matrix_multiply_add(const Matrix& matrix, const SourceView& src, TargetView& tgt, const Configuration& config);

template <typename Matrix, typename SourceView, typename TargetView>
void sparse_matrix_multiply_add(const Matrix& matrix, const SourceView& src, TargetView& tgt, Indexing);

template <typename Matrix, typename SourceView, typename TargetView>
void sparse_matrix_multiply_add(const Matrix& matrix, const SourceView& src, TargetView& tgt, Indexing,
                                const Configuration& config);

class SparseMatrixMultiply {
public:
    SparseMatrixMultiply() = default;
    SparseMatrixMultiply(const std::string& backend): backend_{backend} {}
    SparseMatrixMultiply(const sparse::Backend& backend): backend_(backend) {}

    template <typename Matrix, typename SourceView, typename TargetView>
    void operator()(const Matrix& matrix, const SourceView& src, TargetView& tgt) const {
        multiply(matrix, src, tgt);
    }

    template <typename Matrix, typename SourceView, typename TargetView>
    void operator()(const Matrix& matrix, const SourceView& src, TargetView& tgt, Indexing indexing) const {
        multiply(matrix, src, tgt, indexing);
    }

    template <typename Matrix, typename SourceView, typename TargetView>
    void multiply(const Matrix& matrix, const SourceView& src, TargetView& tgt) const {
        sparse_matrix_multiply(matrix, src, tgt, backend());
    }

    template <typename Matrix, typename SourceView, typename TargetView>
    void multiply(const Matrix& matrix, const SourceView& src, TargetView& tgt, Indexing indexing) const {
        sparse_matrix_multiply(matrix, src, tgt, indexing, backend());
    }

    template <typename Matrix, typename SourceView, typename TargetView>
    void multiply_add(const Matrix& matrix, const SourceView& src, TargetView& tgt) const {
        sparse_matrix_multiply_add(matrix, src, tgt, backend());
    }

    template <typename Matrix, typename SourceView, typename TargetView>
    void multiply_add(const Matrix& matrix, const SourceView& src, TargetView& tgt, Indexing indexing) const {
        sparse_matrix_multiply_add(matrix, src, tgt, indexing, backend());
    }

    const sparse::Backend& backend() const { return backend_; }

private:
    sparse::Backend backend_;
};

namespace sparse {

// Template class which needs (full or partial) specialization for concrete template parameters
template <typename Backend, Indexing, int Rank,  typename MatrixValue, typename SourceValue, typename TargetValue>
struct SparseMatrixMultiply {
    static void multiply(const SparseMatrixView<MatrixValue>&, const View<SourceValue, Rank>&, View<TargetValue, Rank>&,
                         const Configuration&) {
        throw_NotImplemented("SparseMatrixMultiply needs a template specialization with the implementation", Here());
    }
    static void multiply_add(const SparseMatrixView<MatrixValue>&, const View<SourceValue, Rank>&, View<TargetValue, Rank>&,
                             const Configuration&) {
        throw_NotImplemented("SparseMatrixMultiply needs a template specialization with the implementation", Here());
    }
};
}  // namespace sparse

}  // namespace linalg
}  // namespace atlas

#include "SparseMatrixMultiply.tcc"
#include "SparseMatrixMultiply_EckitLinalg.h"
#include "SparseMatrixMultiply_OpenMP.h"
#include "SparseMatrixMultiply_HicSparse.h"
