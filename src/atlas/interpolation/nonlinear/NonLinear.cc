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


namespace {
std::string config_missing_value( const MissingValue::Config& c ) {
    std::string value;
    ATLAS_ASSERT_MSG( c.get( "missing_value_type", value ), "NonLinear: expecting 'missing_value_type'" );
    return value;
}
}  // namespace


NonLinear::NonLinear( const Config& config ) : missingValue_( config_missing_value( config ), config ) {
    ATLAS_ASSERT_MSG( missingValue_, "NonLinear: missingValue setup" );
}


NonLinear::~NonLinear() = default;


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
