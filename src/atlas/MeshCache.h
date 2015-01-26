/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Tiago Quintino
/// @date Mar 2013

#ifndef atlas_MeshCache_h
#define atlas_MeshCache_h

#include "eckit/memory/NonCopyable.h"
#include "eckit/filesystem/LocalPathName.h"

#include "atlas/Mesh.h"

namespace atlas {

//------------------------------------------------------------------------------------------------------

class MeshCache : private eckit::NonCopyable {

public: // methods

    /// @returns true if found cache
    static bool add( const std::string& key, atlas::Mesh& );

    /// @returns true if addition was succcessful
    static bool get( const std::string& key, atlas::Mesh& );

    static eckit::LocalPathName filename(const std::string& key);

};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
