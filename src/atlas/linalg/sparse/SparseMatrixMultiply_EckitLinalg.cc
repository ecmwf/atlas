/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "SparseMatrixMultiply_EckitLinalg.h"

#include "atlas/library/config.h"

#if ATLAS_ECKIT_HAVE_ECKIT_585
#include "eckit/linalg/LinearAlgebraSparse.h"
#else
#include "eckit/linalg/LinearAlgebra.h"
#endif

#include "eckit/linalg/Matrix.h"
#include "eckit/linalg/Vector.h"

#include "atlas/runtime/Exception.h"

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

}  // namespace

void SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 1, const double, double>::apply(
    const SparseMatrix& W, const View<const double, 1>& src, View<double, 1>& tgt, const Configuration& config) {
    ATLAS_ASSERT(src.contiguous());
    ATLAS_ASSERT(tgt.contiguous());
    eckit::linalg::Vector v_src(src.data(), src.size());
    eckit::linalg::Vector v_tgt(tgt.data(), tgt.size());
    eckit_linalg_backend(config).spmv(W.host_matrix(), v_src, v_tgt);
}

void SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 2, const double, double>::apply(
    const SparseMatrix& W, const View<const double, 2>& src, View<double, 2>& tgt, const Configuration& config) {
    ATLAS_ASSERT(src.contiguous());
    ATLAS_ASSERT(tgt.contiguous());
    ATLAS_ASSERT(src.shape(1) >= W.cols());
    ATLAS_ASSERT(tgt.shape(1) >= W.rows());
    eckit::linalg::Matrix m_src(src.data(), src.shape(1), src.shape(0));
    eckit::linalg::Matrix m_tgt(tgt.data(), tgt.shape(1), tgt.shape(0));
    eckit_linalg_backend(config).spmm(W.host_matrix(), m_src, m_tgt);
}

void SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_left, 1, const double, double>::apply(
    const SparseMatrix& W, const View<const double, 1>& src, View<double, 1>& tgt, const Configuration& config) {
    SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 1, const double, double>::apply(W, src, tgt,
                                                                                                        config);
}

}  // namespace sparse
}  // namespace linalg
}  // namespace atlas
