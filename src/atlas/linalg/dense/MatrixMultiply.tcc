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

#include "MatrixMultiply.h"

#include "atlas/linalg/Indexing.h"
#include "atlas/linalg/Introspection.h"
#include "atlas/linalg/View.h"
#include "atlas/linalg/dense/Backend.h"
#include "atlas/runtime/Exception.h"

#if ATLAS_ECKIT_HAVE_ECKIT_585
#include "eckit/linalg/LinearAlgebraDense.h"
#else
#include "eckit/linalg/LinearAlgebra.h"
#endif

namespace atlas {
namespace linalg {


namespace dense {
namespace {
template <typename Backend>
struct MatrixMultiplyHelper {
    template <typename Mat>
    static void apply( const Mat& A, const Mat& B, Mat& C,
                       const eckit::Configuration& config ) {
        MatrixMultiply<Backend>::apply( A,B,C, config );
    }
};

}
}

template <typename Mat>
void matrix_multiply( const Mat& A, const Mat& B, Mat& C, const eckit::Configuration& config ) {
    std::string type = config.getString( "type", dense::current_backend() );
    if ( type == dense::backend::eckit_linalg::type() ) {
        dense::MatrixMultiplyHelper<dense::backend::eckit_linalg>::apply( A, B, C, config );
    }
    else {
        if( type == "openmp" ) {
            type = "generic";
            dense::MatrixMultiplyHelper<dense::backend::eckit_linalg>::apply( A, B, C, util::Config("backend",type) );
        }
#if ATLAS_ECKIT_HAVE_ECKIT_585
        else if( eckit::linalg::LinearAlgebraDense::hasBackend(type) ) {
#else
        else if( eckit::linalg::LinearAlgebra::hasBackend(type) ) {
#endif
            dense::MatrixMultiplyHelper<dense::backend::eckit_linalg>::apply( A, B, C, util::Config("backend",type) );
        }
        else {
            throw_NotImplemented( "matrix_multiply cannot be performed with unsupported backend [" + type + "]",
                                  Here() );
        }
    }
}

template <typename Matrix>
void matrix_multiply( const Matrix& A, const Matrix& B, Matrix& C ) {
    matrix_multiply( A, B, C, dense::Backend() );
}

}  // namespace linalg
}  // namespace atlas
