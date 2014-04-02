// (C) Copyright 1996-2014 ECMWF.

#include "eckit/filesystem/LocalPathName.h"

#include "atlas/Gmsh.hpp"
#include "atlas/grid/MeshCache.h"

using namespace eckit;

namespace atlas {

//------------------------------------------------------------------------------------------------------

MeshCache::MeshCache()
{
}

std::string MeshCache::filename(const std::string &key) const
{
    std::stringstream ss;
    ss << "cache/mesh/" << key << ".gmsh";
    return ss.str();
}

bool MeshCache::add(const std::string &key, Mesh& mesh) const
{
    const std::string fn = filename(key);

    const std::string fn_tmp = fn + ".tmp"; /// @todo this should be more 'random'

    eckit::LocalPathName fpath(fn);

    if( fpath.exists() ) return false;

    fpath.dirName().mkdir();

    eckit::LocalPathName fpath_tmp(fn_tmp);
    fpath_tmp.unlink();

    Log::info() << "inserting mesh in cache (" << fn << ")" << std::endl;

    atlas::Gmsh::write3dsurf(mesh,fn_tmp);

    LocalPathName::rename(fpath_tmp,fpath);

    return true;
}

Mesh* MeshCache::get( const std::string &key ) const
{
    std::string fn = filename(key);

    eckit::LocalPathName fpath(fn);
    if(!fpath.exists())
        return NULL;

    Log::info() << "found mesh in cache (" << fn << ")" << std::endl;

    return atlas::Gmsh::read(fn);
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

