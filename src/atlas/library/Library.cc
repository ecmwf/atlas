/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Baudouin Raoult
/// @author Tiago Quintino
/// @date   August 2016

#include "atlas/library/Library.h"
#include "atlas/library/config.h"
#include "atlas/library/atlas.h"

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

static Library libatlas;

Library::Library() : eckit::system::Library( std::string("atlas") ) {}

const Library& Library::instance() {
    return libatlas;
}

void Library::initialise(int argc, char **argv) {
    return atlas::init(argc,argv);
}

void Library::finalise() {
    return atlas::finalise();
}

const void* Library::addr() const {
    return this;
}

std::string Library::version() const {
    return atlas::version();
}

std::string Library::gitsha1(unsigned int count) const {
    return atlas::git_sha1(count);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas

