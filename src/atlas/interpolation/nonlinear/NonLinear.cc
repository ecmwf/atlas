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

#include "atlas/runtime/Exception.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


NonLinear::NonLinear( const Config& config ) : missingValue_( 0. ) {
    config.get( "missingValue", missingValue_ );
    ATLAS_ASSERT( missingValue_ == missingValue_ );
}


bool NonLinear::missingValue( const double& value ) const {
    return value == missingValue_;
}


NonLinear::~NonLinear() = default;


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
