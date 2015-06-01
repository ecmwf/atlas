/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Tesselation_h
#define atlas_Tesselation_h

#include "atlas/Mesh.h"

namespace atlas {

class Grid;

//----------------------------------------------------------------------------------------------------------------------

class Tesselation {

public:

  /// generate a mesh by triangulating the convex hull of the 3D points
  static void delaunay_triangulation(Mesh& mesh);

  /// generates the cell centres en each cell
  /// @warning only for triangles ATM
  /// @warning this function must be checked for the correct INDEX translation to Fortran
  static void create_cell_centres(Mesh& mesh);

  static void build_mesh(const Grid&, Mesh&);

  /// Ensures or creates the necessary nodal structure in the mesh object
  static void create_mesh_structure(Mesh& mesh, const size_t nb_nodes);

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas

#endif
