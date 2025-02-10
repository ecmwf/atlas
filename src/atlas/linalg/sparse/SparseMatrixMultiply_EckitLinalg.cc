/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "SparseMatrixMultiply_EckitLinalg.h"

#include "atlas/array.h"
#include "atlas/library/config.h"

#if ATLAS_ECKIT_HAVE_ECKIT_585
#include "eckit/linalg/LinearAlgebraSparse.h"
#else
#include "eckit/linalg/LinearAlgebra.h"
#endif

#include "eckit/linalg/Matrix.h"
#include "eckit/linalg/Vector.h"

#include "atlas/runtime/Exception.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"

namespace atlas {
namespace linalg {
namespace sparse {

namespace {

#if ATLAS_ECKIT_HAVE_ECKIT_585
const eckit::linalg::LinearAlgebraSparse& eckit_linalg_backend(const Configuration& config) {
    std::string backend = "default";
    config.get("backend", backend);

    if (backend == "default") {
        return eckit::linalg::LinearAlgebraSparse::backend();
    }
    ATLAS_ASSERT(eckit::linalg::LinearAlgebraSparse::hasBackend(backend));
    return eckit::linalg::LinearAlgebraSparse::getBackend(backend);
}
#else
const eckit::linalg::LinearAlgebra& eckit_linalg_backend(const Configuration& config) {
    std::string backend = "default";
    config.get("backend", backend);
    if (backend == "default") {
        return eckit::linalg::LinearAlgebra::backend();
    }
    ATLAS_ASSERT(eckit::linalg::LinearAlgebra::hasBackend(backend));
    return eckit::linalg::LinearAlgebra::getBackend(backend);
}
#endif

template <typename Value, int Rank>
auto linalg_make_view(atlas::array::ArrayT<double>& array) {
    auto v_array = array::make_view<Value, Rank>(array);
    return atlas::linalg::make_view(v_array);
}

}  // namespace

void SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 1, double, const double, double>::multiply(
    const SparseMatrixView<double>& W, const View<const double, 1>& src, View<double, 1>& tgt, const Configuration& config) {
    ATLAS_ASSERT(src.contiguous());
    ATLAS_ASSERT(tgt.contiguous());
    eckit::linalg::Vector v_src(src.data(), src.size());
    eckit::linalg::Vector v_tgt(tgt.data(), tgt.size());
    eckit::linalg::SparseMatrix m_W = make_non_owning_eckit_sparse_matrix(W);
    eckit_linalg_backend(config).spmv(m_W, v_src, v_tgt);
}

void SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 2, double, const double, double>::multiply(
    const SparseMatrixView<double>& W, const View<const double, 2>& src, View<double, 2>& tgt, const Configuration& config) {
    ATLAS_ASSERT(src.contiguous());
    ATLAS_ASSERT(tgt.contiguous());
    ATLAS_ASSERT(src.shape(1) >= W.cols());
    ATLAS_ASSERT(tgt.shape(1) >= W.rows());
    eckit::linalg::Matrix m_src(src.data(), src.shape(1), src.shape(0));
    eckit::linalg::Matrix m_tgt(tgt.data(), tgt.shape(1), tgt.shape(0));
    eckit::linalg::SparseMatrix m_W = make_non_owning_eckit_sparse_matrix(W);
    eckit_linalg_backend(config).spmm(m_W, m_src, m_tgt);
}

void SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_left, 1, double, const double, double>::multiply(
    const SparseMatrixView<double>& W, const View<const double, 1>& src, View<double, 1>& tgt, const Configuration& config) {
    SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 1, double, const double, double>::multiply(W, src, tgt,
                                                                                                        config);
}

void SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 1, double, const double, double>::multiply_add(
    const SparseMatrixView<double>& W, const View<const double, 1>& src, View<double, 1>& tgt, const Configuration& config) {

    array::ArrayT<double> tmp(src.shape(0));
    auto v_tmp = linalg_make_view<double, 1>(tmp);
    v_tmp.assign(0.);

    SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 1, double, const double, double>::multiply(W, src, v_tmp, config);

    for (idx_t t = 0; t < tmp.shape(0); ++t) {
        tgt(t) += v_tmp(t);
    }
}

void SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 2, double, const double, double>::multiply_add(
    const SparseMatrixView<double>& W, const View<const double, 2>& src, View<double, 2>& tgt, const Configuration& config) {

    array::ArrayT<double> tmp(src.shape(0), src.shape(1));
    auto v_tmp = linalg_make_view<double, 2>(tmp);
    v_tmp.assign(0.);

    SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 2, double, const double, double>::multiply(W, src, v_tmp, config);

    for (idx_t t = 0; t < tmp.shape(0); ++t) {
        for (idx_t k = 0; k < tmp.shape(1); ++k) {
            tgt(t, k) += v_tmp(t, k);
        }
    }
}

void SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_left, 1, double, const double, double>::multiply_add(
    const SparseMatrixView<double>& W, const View<const double, 1>& src, View<double, 1>& tgt, const Configuration& config) {
    SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 1, double, const double, double>::multiply_add(W, src, tgt, config);
}

}  // namespace sparse
}  // namespace linalg
}  // namespace atlas
