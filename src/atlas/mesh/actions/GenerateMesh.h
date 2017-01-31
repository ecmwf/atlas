/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_GenerateMesh_h
#define atlas_GenerateMesh_h

#include "atlas/grid/Structured.h"
#include "atlas/mesh/Mesh.h"

namespace atlas {
namespace grid {
class Structured;
class GridDistribution;
}
namespace mesh {
class Mesh;
}
}

namespace atlas {
namespace mesh {
namespace actions {

//mesh::Mesh* generate_mesh (const grid::Structured& grid);

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

#define STRUCTUREDGRID grid::Structured
#define grid_GridDistribution grid::GridDistribution
#define mesh_Mesh mesh::Mesh
extern "C"
{
  mesh_Mesh* atlas__generate_mesh (STRUCTUREDGRID* grid);
  mesh_Mesh* atlas__generate_mesh_with_distribution (STRUCTUREDGRID* grid, grid_GridDistribution* distribution);
}

#undef grid_GridDistribution
#undef STRUCTUREDGRID

// ------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas

#endif
