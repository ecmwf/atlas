/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_GenerateMesh_h
#define atlas_GenerateMesh_h

#include "atlas/Mesh.h"
#include "atlas/grids/ReducedGrid.h"

namespace atlas {

class GridDistribution;
class Mesh;

namespace grids {
  class ReducedGrid;
}

namespace actions {

Mesh* generate_mesh (const grids::ReducedGrid& rgg);

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

typedef grids::ReducedGrid __ReducedGrid;
extern "C"
{
  Mesh* atlas__generate_mesh (__ReducedGrid* rgg);
  Mesh* atlas__generate_mesh_with_distribution (__ReducedGrid* rgg, GridDistribution* distribution);

}

// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas

#endif
