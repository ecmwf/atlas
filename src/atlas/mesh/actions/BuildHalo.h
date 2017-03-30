/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef BuildHalo_h
#define BuildHalo_h

#include <string>

namespace atlas {
  class Mesh;

namespace mesh {
namespace actions {

/// @brief Enlarge each partition of the mesh with a halo of elements
/// @param [inout] mesh      The mesh to enlarge
/// @param [in]    nb_elems  Size of the halo
/// @author Willem Deconinck
/// @date June 2014
void build_halo( Mesh& mesh, int nb_elems );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  void atlas__build_halo (Mesh::mesh_t* mesh, int nb_elems);
}
// ------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas

#endif // BuildHalo_h
