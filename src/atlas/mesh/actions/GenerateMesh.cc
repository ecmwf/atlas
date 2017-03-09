/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grid/Distribution.h"
#include "atlas/mesh/actions/GenerateMesh.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/grid/Grid.h"

namespace atlas {
namespace mesh {
namespace actions {

// ------------------------------------------------------------------

namespace {
Mesh* generate_mesh(const grid::Grid& rgg)
{
  Log::info() << "Deprecated function [generate_mesh] used.\n"
              << "Consider using meshgenerator::Structured directly."
              << std::endl;

  meshgenerator::StructuredMeshGenerator generate;
  return generate(rgg);
}
}

// ------------------------------------------------------------------


Mesh* atlas__generate_mesh(grid::Grid::grid_t* rgg)
{
  ATLAS_ERROR_HANDLING( return generate_mesh( grid::Grid(rgg) ); );
  return NULL;
}

// ------------------------------------------------------------------


Mesh* atlas__generate_mesh_with_distribution(grid::Grid::grid_t* rgg, grid::Distribution::impl_t* distribution)
{
  ATLAS_ERROR_HANDLING(
        meshgenerator::StructuredMeshGenerator generate;
        return generate( grid::Grid(rgg), grid::Distribution(distribution) );
  );
  return NULL;
}

// ------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas
