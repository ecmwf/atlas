/*
 * (C) Copyright 1996-2016 ECMWF.
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

#include "eckit/system/Library.h"

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

class LibAtlas : public eckit::system::Library {
public:
    LibAtlas() : Library("atlas") {}
protected:
    const void* addr() const { return this; }
};

static LibAtlas libatlas;

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas

