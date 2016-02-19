/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grid/GridDistribution.h"
#include "atlas/mesh/actions/GenerateMesh.h"
#include "atlas/mesh/generators/ReducedGridMeshGenerator.h"
#include "atlas/util/runtime/Log.h"
#include "atlas/util/runtime/ErrorHandling.h"
#include "atlas/util/parallel/mpi/mpi.h"

using namespace atlas::meshgen;
using namespace atlas::grids;
namespace atlas {
namespace actions {

// ------------------------------------------------------------------

Mesh* generate_mesh (const ReducedGrid& rgg)
{
  Log::info() << "Deprecated function [generate_mesh] used.\n"
              << "Consider using ReducedGridMeshGenerator directly."
              << std::endl;

  ReducedGridMeshGenerator generate;
  generate.options.set( "nb_parts", eckit::mpi::size() );
  generate.options.set( "part"    , eckit::mpi::rank() );
  return generate(rgg);
}

// ------------------------------------------------------------------


Mesh* atlas__generate_mesh (ReducedGrid* rgg)
{
  ATLAS_ERROR_HANDLING( return generate_mesh(*rgg); );
  return NULL;
}

// ------------------------------------------------------------------


Mesh* atlas__generate_mesh_with_distribution (ReducedGrid* rgg, GridDistribution* distribution)
{
  ATLAS_ERROR_HANDLING(
        ReducedGridMeshGenerator generate;
        return generate(*rgg, *distribution);
  );
  return NULL;
}

// ------------------------------------------------------------------


} // namespace actions
} // namespace atlas
