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


#include "atlas/interpolation/nonlinear/NonLinear.h"

#include <cmath>

#include "atlas/runtime/Exception.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


/// @brief CompareNaN Indicate missing value if NaN
struct CompareNaN : NonLinear::Compare {
    bool operator()( const double& value ) const override { return std::isnan( value ); }
};


/// @brief CompareValue Indicate missing value if compares equally to pre-defined value
struct CompareValue : NonLinear::Compare {
    CompareValue( double missingValue ) : missingValue_( missingValue ) { ATLAS_ASSERT( !std::isnan( missingValue ) ); }

    bool operator()( const double& value ) const override { return value == missingValue_; }

    double missingValue_;
};


NonLinear::NonLinear( const Config& config ) : missingValue_( new CompareNaN() ) {
    double missingValue;
    if ( config.get( "missingValue", missingValue ) ) {
        missingValue_.reset( new CompareValue( missingValue ) );
    }
}


NonLinear::~NonLinear() = default;


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
