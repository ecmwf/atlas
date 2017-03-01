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

#include "atlas/runtime/LibAtlas.h"

#include "atlas/internals/atlas_version.h"

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

LibAtlas::LibAtlas() : Library( std::string("atlas") ) {}

const LibAtlas& LibAtlas::instance()
{
    static const LibAtlas libatlas;
    return libatlas;
}

const void* LibAtlas::addr() const { return this; }

std::string LibAtlas::version() const { return atlas_version_str(); }

std::string LibAtlas::gitsha1(unsigned int count) const {
    return atlas_git_sha1_abbrev(count);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas

