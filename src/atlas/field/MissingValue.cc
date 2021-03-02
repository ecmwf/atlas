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


#include "atlas/field/MissingValue.h"

#include <cmath>

#include "eckit/types/FloatCompare.h"

#include "atlas/field/Field.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Metadata.h"


namespace atlas {
namespace field {


namespace {
std::string config_type( const MissingValue::Config& c ) {
    std::string value;
    c.get( "type", value );
    return value;
}


std::string field_type( const Field& f ) {
    std::string value;
    if ( f.metadata().get( "missing_value_type", value ) ) {
        value += "-" + f.datatype().str();
    }
    return value;
}
}  // namespace


MissingValue::MissingValue() : Handle( nullptr ) {}


MissingValue::MissingValue( const MissingValue::Config& config ) :
    Handle( detail::MissingValueFactory::build( config_type( config ), config ) ) {}


MissingValue::MissingValue( const std::string& type, const MissingValue::Config& config ) :
    Handle( detail::MissingValueFactory::build( type, config ) ) {}


MissingValue::MissingValue( const Field& field ) :
    Handle( detail::MissingValueFactory::build( field_type( field ), field.metadata() ) ) {}


MissingValue::MissingValue( const std::string& type, const Field& field ) :
    Handle( detail::MissingValueFactory::build( type, field.metadata() ) ) {}


bool MissingValue::operator()( const double& value ) const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue::operator()( const double& ) ObjectHandle not setup" );
    return get()->operator()( value );
}


bool MissingValue::operator()( const float& value ) const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue::operator()( const float& ): ObjectHandle not setup" );
    return get()->operator()( value );
}


bool MissingValue::operator()( const int& value ) const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue::operator()( const int& ): ObjectHandle not setup" );
    return get()->operator()( value );
}


bool MissingValue::operator()( const long& value ) const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue::operator()( const long& ): ObjectHandle not setup" );
    return get()->operator()( value );
}


bool MissingValue::operator()( const unsigned long& value ) const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue::operator()( const unsigned long& ): ObjectHandle not setup" );
    return get()->operator()( value );
}


bool MissingValue::isnan() const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue: ObjectHandle not setup" );
    return get()->isnan();
}


const detail::MissingValue& MissingValue::ref() const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue: ObjectHandle not setup" );
    return *get();  // (one dereferencing level less)
}


void MissingValue::metadata( Field& field ) const {
    ATLAS_ASSERT_MSG( operator bool(), "MissingValue: ObjectHandle not setup" );
    get()->metadata( field );
}


}  // namespace field
}  // namespace atlas
