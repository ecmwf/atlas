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
#include "eckit/linalg/Matrix.h"

#include "atlas/linalg/Indexing.h"
#include "atlas/linalg/View.h"
#include "atlas/linalg/dense/Backend.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace linalg {

using Matrix        = eckit::linalg::Matrix;
using Configuration = eckit::Configuration;

// C = A . B
template <typename Matrix>
void matrix_multiply(const Matrix& A, const Matrix& B, Matrix& C);

template <typename Matrix>
void matrix_multiply(const Matrix& A, const Matrix& B, Matrix& C, const eckit::Configuration&);

class MatrixMultiply {
public:
    MatrixMultiply() = default;
    MatrixMultiply(const std::string& backend): backend_{backend} {}
    MatrixMultiply(const dense::Backend& backend): backend_(backend) {}

    void operator()(const Matrix& A, const Matrix& B, Matrix& C) const { matrix_multiply(A, B, C, backend()); }

    const dense::Backend& backend() const { return backend_; }

private:
    dense::Backend backend_;
};

namespace dense {

// Template class which needs (full or partial) specialization for concrete template parameters
template <typename Backend>
struct MatrixMultiply {
    static void apply(const Matrix& A, const Matrix& B, Matrix& C, const Configuration&) {
        throw_NotImplemented("MatrixMultiply needs a template specialization with the implementation", Here());
    }
};
}  // namespace dense

}  // namespace linalg
}  // namespace atlas

#include "MatrixMultiply.tcc"
#include "MatrixMultiply_EckitLinalg.h"
