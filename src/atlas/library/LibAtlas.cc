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

#include "atlas/library/LibAtlas.h"
#include "atlas/library/config.h"

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

static LibAtlas libatlas;

LibAtlas::LibAtlas() : Library( std::string("atlas") ) {}

const LibAtlas& LibAtlas::instance() {
    return libatlas;
}

const void* LibAtlas::addr() const {
    return this;
}

std::string LibAtlas::version() const {
    return atlas::version();
}

std::string LibAtlas::gitsha1(unsigned int count) const {
    return atlas::git_sha1(count);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas

