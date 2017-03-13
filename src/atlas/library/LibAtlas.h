/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "eckit/system/Library.h"

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

class LibAtlas : public eckit::system::Library {
public:

    LibAtlas();

    static const LibAtlas& instance();

protected:

    const void* addr() const;

    virtual std::string version() const;

    virtual std::string gitsha1(unsigned int count) const;

};

typedef LibAtlas Atlas;

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas

