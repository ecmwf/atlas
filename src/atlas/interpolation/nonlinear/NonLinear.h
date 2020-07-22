/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include "eckit/config/Parametrisation.h"
#include "eckit/linalg/SparseMatrix.h"


namespace atlas {
class Field;
}


namespace atlas {
namespace interpolation {
namespace nonlinear {


/**
 * @brief NonLinear class applies non-linear corrections to an interpolation matrix, given a field with missing values.
 * The interpolatation are re-weighted to factor those values out of the resulting field.
 */
class NonLinear {
public:
    using Config = eckit::Parametrisation;
    using Matrix = eckit::linalg::SparseMatrix;
    using Scalar = eckit::linalg::Scalar;
    using Size   = eckit::linalg::Size;

    /**
     * @brief NonLinear ctor
     * @param [in] config with optional "missingValue", default NaN indicates missing value in field
     */
    NonLinear( const Config& config );

    /**
     * @brief NonLinear dtor
     */
    virtual ~NonLinear();

    /**
     * @brief Apply non-linear corrections to interpolation matrix
     * @param [inout] W interpolation matrix
     * @param [in] f field with missing values information
     * @return if W was modified
     */
    virtual bool execute( Matrix& W, const Field& f ) const = 0;

    /**
     * @brief Return if value is considered a missing value
     * @param [in] value to evaluate
     * @return if value is considered a missing value
     */
    bool missingValue( const double& ) const;

private:
    double missingValue_;
    const bool missingValueIsNaN_;
};


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
