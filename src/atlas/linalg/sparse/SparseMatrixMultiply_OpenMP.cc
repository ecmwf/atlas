/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/linalg/sparse/SparseMatrixMultiply_OpenMP.h"

#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace linalg {
namespace sparse {

template <bool SetZero, typename SourceValue, typename TargetValue>
void spmv_layout_left(const SparseMatrix& W, const View<SourceValue, 1>& src, View<TargetValue, 1>& tgt) {
    using Value       = TargetValue;
    const auto outer  = W.outer();
    const auto index  = W.inner();
    const auto weight = W.data();
    const idx_t rows  = static_cast<idx_t>(W.rows());

    ATLAS_ASSERT(src.shape(0) >= W.cols());
    ATLAS_ASSERT(tgt.shape(0) >= W.rows());

    atlas_omp_parallel_for(idx_t r = 0; r < rows; ++r) {
        if constexpr (SetZero) {
            tgt[r] = 0.;
        }
        for (idx_t c = outer[r]; c < outer[r + 1]; ++c) {
            idx_t n = index[c];
            Value w = static_cast<Value>(weight[c]);
            tgt[r] += w * src[n];
        }
    }
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 1, SourceValue, TargetValue>::multiply(
    const SparseMatrix& W, const View<SourceValue, 1>& src, View<TargetValue, 1>& tgt, const Configuration&) {
    spmv_layout_left<true>(W, src, tgt);
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 1, SourceValue, TargetValue>::multiply_add(
    const SparseMatrix& W, const View<SourceValue, 1>& src, View<TargetValue, 1>& tgt, const Configuration&) {
    spmv_layout_left<false>(W, src, tgt);
}

template <bool SetZero, typename SourceValue, typename TargetValue>
void spmm_layout_left(const SparseMatrix& W, const View<SourceValue, 2>& src, View<TargetValue, 2>& tgt) {
    using Value       = TargetValue;
    const auto outer  = W.outer();
    const auto index  = W.inner();
    const auto weight = W.data();
    const idx_t rows  = static_cast<idx_t>(W.rows());
    const idx_t Nk    = src.shape(1);

    ATLAS_ASSERT(src.shape(0) >= W.cols());
    ATLAS_ASSERT(tgt.shape(0) >= W.rows());

    atlas_omp_parallel_for(idx_t r = 0; r < rows; ++r) {
        if constexpr (SetZero) {
            for (idx_t k = 0; k < Nk; ++k) {
                tgt(r, k) = 0.;
            }
        }
        for (idx_t c = outer[r]; c < outer[r + 1]; ++c) {
            idx_t n = index[c];
            Value w = static_cast<Value>(weight[c]);
            for (idx_t k = 0; k < Nk; ++k) {
                tgt(r, k) += w * src(n, k);
            }
        }
    }
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 2, SourceValue, TargetValue>::multiply(
    const SparseMatrix& W, const View<SourceValue, 2>& src, View<TargetValue, 2>& tgt, const Configuration&) {
    spmm_layout_left<true>(W, src, tgt);
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 2, SourceValue, TargetValue>::multiply_add(
    const SparseMatrix& W, const View<SourceValue, 2>& src, View<TargetValue, 2>& tgt, const Configuration&) {
    spmm_layout_left<false>(W, src, tgt);
}

template <bool SetZero, typename SourceValue, typename TargetValue>
void spmt_layout_left(const SparseMatrix& W, const View<SourceValue, 3>& src, View<TargetValue, 3>& tgt) {
    if (src.contiguous() && tgt.contiguous()) {
        // We can take a more optimized route by reducing rank
        auto src_v = View<SourceValue, 2>(src.data(), array::make_shape(src.shape(0), src.stride(0)));
        auto tgt_v = View<TargetValue, 2>(tgt.data(), array::make_shape(tgt.shape(0), tgt.stride(0)));
        spmm_layout_left<SetZero>(W, src_v, tgt_v);
        return;
    }
    using Value       = TargetValue;
    const auto outer  = W.outer();
    const auto index  = W.inner();
    const auto weight = W.data();
    const idx_t rows  = static_cast<idx_t>(W.rows());
    const idx_t Nk    = src.shape(1);
    const idx_t Nl    = src.shape(2);

    atlas_omp_parallel_for(idx_t r = 0; r < rows; ++r) {
        if constexpr (SetZero) {
            for (idx_t k = 0; k < Nk; ++k) {
                for (idx_t l = 0; l < Nl; ++l) {
                    tgt(r, k, l) = 0.;
                }
            }
        }
        for (idx_t c = outer[r]; c < outer[r + 1]; ++c) {
            idx_t n       = index[c];
            const Value w = static_cast<Value>(weight[c]);
            for (idx_t k = 0; k < Nk; ++k) {
                for (idx_t l = 0; l < Nl; ++l) {
                    tgt(r, k, l) += w * src(n, k, l);
                }
            }
        }
    }
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 3, SourceValue, TargetValue>::multiply(
    const SparseMatrix& W, const View<SourceValue, 3>& src, View<TargetValue, 3>& tgt, const Configuration& config) {
    spmt_layout_left<true>(W, src, tgt);
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 3, SourceValue, TargetValue>::multiply_add(
    const SparseMatrix& W, const View<SourceValue, 3>& src, View<TargetValue, 3>& tgt, const Configuration& config) {
    spmt_layout_left<false>(W, src, tgt);
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 1, SourceValue, TargetValue>::multiply(
    const SparseMatrix& W, const View<SourceValue, 1>& src, View<TargetValue, 1>& tgt, const Configuration& config) {
    return SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 1, SourceValue, TargetValue>::multiply(W, src, tgt, config);
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 1, SourceValue, TargetValue>::multiply_add(
    const SparseMatrix& W, const View<SourceValue, 1>& src, View<TargetValue, 1>& tgt, const Configuration& config) {
    return SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 1, SourceValue, TargetValue>::multiply_add(W, src, tgt, config);
}

template <bool SetZero, typename SourceValue, typename TargetValue>
void spmm_layout_right(const SparseMatrix& W, const View<SourceValue, 2>& src, View<TargetValue, 2>& tgt) {
    using Value       = TargetValue;
    const auto outer  = W.outer();
    const auto index  = W.inner();
    const auto weight = W.data();
    const idx_t rows  = static_cast<idx_t>(W.rows());
    const idx_t Nk    = src.shape(0);

    ATLAS_ASSERT(src.shape(1) >= W.cols());
    ATLAS_ASSERT(tgt.shape(1) >= W.rows());

    atlas_omp_parallel_for(idx_t r = 0; r < rows; ++r) {
        if constexpr (SetZero) {
            for (idx_t k = 0; k < Nk; ++k) {
                tgt(k, r) = 0.;
            }
        }
        for (idx_t c = outer[r]; c < outer[r + 1]; ++c) {
            idx_t n = index[c];
            Value w = static_cast<Value>(weight[c]);
            for (idx_t k = 0; k < Nk; ++k) {
                tgt(k, r) += w * src(k, n);
            }
        }
    }
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 2, SourceValue, TargetValue>::multiply(
    const SparseMatrix& W, const View<SourceValue, 2>& src, View<TargetValue, 2>& tgt, const Configuration&) {
    spmm_layout_right<true>(W, src, tgt);
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 2, SourceValue, TargetValue>::multiply_add(
    const SparseMatrix& W, const View<SourceValue, 2>& src, View<TargetValue, 2>& tgt, const Configuration&) {
    spmm_layout_right<false>(W, src, tgt);
}

template <bool SetZero, typename SourceValue, typename TargetValue>
void spmt_layout_right(const SparseMatrix& W, const View<SourceValue, 3>& src, View<TargetValue, 3>& tgt) {
    if (src.contiguous() && tgt.contiguous()) {
        // We can take a more optimized route by reducing rank
        auto src_v = View<SourceValue, 2>(src.data(), array::make_shape(src.shape(0), src.stride(0)));
        auto tgt_v = View<TargetValue, 2>(tgt.data(), array::make_shape(tgt.shape(0), tgt.stride(0)));
        spmm_layout_right<SetZero>(W, src_v, tgt_v);
        return;
    }
    using Value       = TargetValue;
    const auto outer  = W.outer();
    const auto index  = W.inner();
    const auto weight = W.data();
    const idx_t rows  = static_cast<idx_t>(W.rows());
    const idx_t Nk    = src.shape(1);
    const idx_t Nl    = src.shape(0);

    atlas_omp_parallel_for(idx_t r = 0; r < rows; ++r) {
        if constexpr (SetZero){
            for (idx_t k = 0; k < Nk; ++k) {
                for (idx_t l = 0; l < Nl; ++l) {
                    tgt(l, k, r) = 0.;
                }
            }
        }
        for (idx_t c = outer[r]; c < outer[r + 1]; ++c) {
            idx_t n       = index[c];
            const Value w = static_cast<Value>(weight[c]);
            for (idx_t k = 0; k < Nk; ++k) {
                for (idx_t l = 0; l < Nl; ++l) {
                    tgt(l, k, r) += w * src(l, k, n);
                }
            }
        }
    }
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 3, SourceValue, TargetValue>::multiply(
    const SparseMatrix& W, const View<SourceValue, 3>& src, View<TargetValue, 3>& tgt, const Configuration& config) {
    spmt_layout_right<true>(W, src, tgt);
}

template <typename SourceValue, typename TargetValue>
void SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 3, SourceValue, TargetValue>::multiply_add(
    const SparseMatrix& W, const View<SourceValue, 3>& src, View<TargetValue, 3>& tgt, const Configuration& config) {
    spmt_layout_right<false>(W, src, tgt);
}

#define EXPLICIT_TEMPLATE_INSTANTIATION(TYPE)                                                           \
    template struct SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 1, TYPE const, TYPE>;  \
    template struct SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 2, TYPE const, TYPE>;  \
    template struct SparseMatrixMultiply<backend::openmp, Indexing::layout_left, 3, TYPE const, TYPE>;  \
    template struct SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 1, TYPE const, TYPE>; \
    template struct SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 2, TYPE const, TYPE>; \
    template struct SparseMatrixMultiply<backend::openmp, Indexing::layout_right, 3, TYPE const, TYPE>;

EXPLICIT_TEMPLATE_INSTANTIATION(double);
EXPLICIT_TEMPLATE_INSTANTIATION(float);

}  // namespace sparse
}  // namespace linalg
}  // namespace atlas
