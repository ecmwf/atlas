/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <string>

#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/config/Resource.h"

#include "atlas/atlas_version.h"
#include "atlas/io/Gmsh.h"
#include "atlas/MeshCache.h"

using eckit::Log;
using eckit::PathName;
using eckit::Resource;
using atlas::io::Gmsh;

namespace atlas {

//------------------------------------------------------------------------------------------------------

MeshCache::MeshCache() : CacheManager("atlas/mesh") {}

const char* MeshCache::version() const { return atlas_version_str(); }
const char* MeshCache::extension() const { return ".msh"; }

std::string MeshCache::generate_key(const Grid& g) const {
  std::ostringstream s;
  s << g.unique_id();
  return s.str();
}

void MeshCache::print(std::ostream &s) const
{
  s << "MeshCache[";
  CacheManager::print(s);
  s << "name=" << name() << ","
    << "version=" << version() << ","
    << "extention=" << extension() << ","
    << "]";
}

void MeshCache::insert(const Grid& grid, const Mesh& mesh) const {
  key_t key = generate_key(grid);

  PathName tmp_path = stage(key);

  Log::info() << "Inserting mesh in cache (" << tmp_path << ")" << std::endl;

  {
    //      FileHandle f(tmp_path, true);
    //      f.openForWrite(0); AutoClose closer(f);

    /// @TODO : change Gmsh writer to use FileHandle

    Gmsh gmsh;
    gmsh.options.set<std::string>("nodes","xyz");
    gmsh.write(mesh, tmp_path.asString());
  }

  commit(key, tmp_path);
}

bool MeshCache::retrieve(const Grid& grid, Mesh& mesh) const {
  key_t key = generate_key(grid);

  PathName path;

  Log::info() << "Looking for cached mesh for grid (" << grid.unique_id() << ")" << std::endl;

  if (!get(key, path)) return false;

  Log::info() << "Found mesh in cache (" << path << ")" << std::endl;

  {
    //      FileHandle f(path, true);
    //      f.openForWrite(0); AutoClose closer(f);

    /// @TODO : change Gmsh reader to use FileHandle

    Gmsh gmsh;
    gmsh.read(path, mesh);
  }

  return true;
}

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
