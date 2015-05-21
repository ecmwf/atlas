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

#include <string>

#include "eckit/container/CacheManager.h"
#include "eckit/filesystem/PathName.h"

#include "atlas/Mesh.h"
#include "atlas/Grid.h"

namespace atlas {

//------------------------------------------------------------------------------------------------------

class MeshCache  : public eckit::CacheManager {

 public: // methods

  explicit MeshCache(bool nocache = false);

  /// Tries to retrieve a cached Mesh
  /// @returns true if found cache
  bool retrieve(const atlas::Grid& g, atlas::Mesh& m) const;

  /// Inserts a cached WeightMatrix, overwritting any existing entry
  /// @returns true if insertion successful cache
  void insert(const atlas::Grid& g, const atlas::Mesh& m) const;

 protected:

  virtual void print(std::ostream& s) const;

 private:

  const bool nocache_;

  std::string generate_key(const atlas::Grid& g) const;

  virtual const char* version() const;
  virtual const char* extension() const;

};

//------------------------------------------------------------------------------------------------------

}  // namespace atlas

#endif
