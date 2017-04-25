/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_mesh_Halo_h
#define atlas_mesh_Halo_h

#include <cstddef>

namespace atlas {
    class Mesh;
}

namespace atlas {
namespace mesh {

// -------------------------------------------------------------------

class Halo
{
public:
  Halo(const Mesh& mesh);
  Halo(const size_t size) : size_(size) {}
  size_t size() const { return size_; }
private:
  size_t size_;
};

} // namespace mesh
} // namespace atlas

#endif // atlas_mesh_Halo_h
