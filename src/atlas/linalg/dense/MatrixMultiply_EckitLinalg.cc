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

#include "MatrixMultiply_EckitLinalg.h"

#include "eckit/linalg/LinearAlgebra.h"
#include "eckit/linalg/Matrix.h"
#include "eckit/linalg/Vector.h"

#include "atlas/runtime/Exception.h"

namespace atlas {
namespace linalg {
namespace dense {

namespace {
const eckit::linalg::LinearAlgebra& eckit_linalg_backend(const Configuration& config) {
    std::string backend = "default";
    config.get("backend", backend);
    if (backend == "default") {
        return eckit::linalg::LinearAlgebra::backend();
    }
    ATLAS_ASSERT(eckit::linalg::LinearAlgebra::hasBackend(backend));
    return eckit::linalg::LinearAlgebra::getBackend(backend);
}

}  // namespace

void MatrixMultiply<backend::eckit_linalg>::apply(const Matrix& A, const Matrix& B, Matrix& C,
                                                  const Configuration& config) {
    eckit_linalg_backend(config).gemm(A, B, C);
}


}  // namespace dense
}  // namespace linalg
}  // namespace atlas
