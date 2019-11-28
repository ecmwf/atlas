/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string>

#include "atlas/output/detail/GmshImpl.h"
#include "atlas/output/detail/GmshInterface.h"

namespace atlas {
namespace output {
namespace detail {

// -----------------------------------------------------------------------------

extern "C" {

GmshImpl* atlas__output__Gmsh__create_pathname_mode( const char* pathname, const char* mode ) {
    return new GmshImpl( std::string( pathname ), std::string( mode ) );
}
GmshImpl* atlas__output__Gmsh__create_pathname_mode_config( const char* pathname, const char* mode,
                                                            const eckit::Parametrisation* config ) {
    return new GmshImpl( std::string( pathname ), std::string( mode ), *config );
}

}  // extern C

//----------------------------------------------------------------------------------------------------------------------
}  // namespace detail
}  // namespace output
}  // namespace atlas
