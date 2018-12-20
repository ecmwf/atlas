/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/exception/Exceptions.h"

#include "atlas/runtime/Exception.h"

namespace atlas {

void throw_NotImplemented( const eckit::CodeLocation& loc ) {
    throw eckit::NotImplemented( loc );
}

void throw_NotImplemented( const std::string& msg, const eckit::CodeLocation& loc ) {
    throw eckit::NotImplemented( msg, loc );
}

void throw_AssertionFailed( const std::string& msg ) {
    throw eckit::AssertionFailed( msg );
}

void throw_AssertionFailed( const std::string& msg, const eckit::CodeLocation& loc ) {
    throw eckit::AssertionFailed( msg, loc );
}

void throw_Exception( const std::string& msg ) {
    throw eckit::Exception( msg );
}

void throw_Exception( const std::string& msg, const eckit::CodeLocation& loc ) {
    throw eckit::Exception( msg, loc );
}

void throw_SeriousBug( const std::string& msg ) {
    throw eckit::SeriousBug( msg );
}

void throw_SeriousBug( const std::string& msg, const eckit::CodeLocation& loc ) {
    throw eckit::SeriousBug( msg, loc );
}

}  // namespace atlas
