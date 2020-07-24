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


#include "atlas/interpolation/MissingValue.h"

#include <cmath>

#include "eckit/types/FloatCompare.h"

#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"


namespace atlas {
namespace interpolation {


namespace {
using Config = MissingValue::Config;


std::string config_type( const Config& c ) {
    std::string value;
    ATLAS_ASSERT( c.get( "type", value ) );
    return value;
}
}  // namespace


MissingValue::MissingValue() : Handle( nullptr ) {}


MissingValue::MissingValue( const MissingValue::Config& config ) :
    Handle( nonlinear::MissingValueFactory::build( config_type( config ), config ) ) {}


MissingValue::MissingValue( const std::string& type, const MissingValue::Config& config ) :
    Handle( nonlinear::MissingValueFactory::build( type, config ) ) {}


bool MissingValue::operator()( const double& value ) const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue: ObjectHandle not setup" );
    return get()->operator()( value );
}


bool MissingValue::isnan() const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue: ObjectHandle not setup" );
    return get()->isnan();
}


const nonlinear::MissingValue& MissingValue::ref() const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue: ObjectHandle not setup" );
    return *get();  // (one dereferencing level less)
}


}  // namespace interpolation
}  // namespace atlas
