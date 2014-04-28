// (C) Copyright 1996-2014 ECMWF.

#include "eckit/filesystem/LocalPathName.h"

#include "atlas/Gmsh.hpp"
#include "atlas/grid/MeshCache.h"

using namespace eckit;

namespace atlas {

//------------------------------------------------------------------------------------------------------

std::string MeshCache::filename(const std::string &key)
{
    std::stringstream ss;
    ss << "cache/atlas/mesh/" << key << ".gmsh";
    return ss.str();
}

bool MeshCache::add(const std::string& key, Mesh& mesh)
{
    LocalPathName file( filename(key) );

    if( file.exists() )
    {
        Log::debug() << "MeshCache entry " << file << " already exists ..." << std::endl;
        return false;
    }

    file.dirName().mkdir();  // ensure directory exists

    // unique file name avoids race conditions on the file from multiple processes

    LocalPathName tmpfile ( LocalPathName::unique(file) );

    Log::info() << "inserting mesh in cache (" << file << ")" << std::endl;

    atlas::Gmsh::write3dsurf(mesh,tmpfile);

    // now try to rename the file to its file pathname

    try
    {
        LocalPathName::rename( tmpfile, file );
    }
    catch( FailedSystemCall& e ) // ignore failed system call -- another process nay have created the file meanwhile
    {
        Log::debug() << "Failed rename of cache file -- " << e.what() << std::endl;
    }

    return true;
}

bool MeshCache::get(const std::string &key, Mesh& mesh)
{
    LocalPathName file( filename(key) );

    if( ! file.exists() )
    {
        return false;
    }

    Log::info() << "found mesh in cache (" << file << ")" << std::endl;

    atlas::Gmsh::read(file,mesh);

    return true;
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

