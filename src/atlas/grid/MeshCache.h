// (C) Copyright 1996-2014 ECMWF.

#ifndef atlas_grid_MeshCache_hpp
#define atlas_grid_MeshCache_hpp

#include "eckit/memory/NonCopyable.h"

#include "atlas/Mesh.hpp"

namespace atlas {

//------------------------------------------------------------------------------------------------------

class MeshCache : private eckit::NonCopyable {

public: // methods

    MeshCache();

    bool add( const std::string& key, atlas::Mesh& ) const;
    atlas::Mesh* get( const std::string& key ) const;

    std::string filename(const std::string& key) const;

};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
