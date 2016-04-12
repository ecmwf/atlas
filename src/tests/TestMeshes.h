/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/internals/atlas_config.h"
#include "atlas/grid/grids.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/grid/global/gaussian/ReducedGaussian.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/parallel/mpi/mpi.h"

using namespace atlas;
using namespace atlas::grid;
using namespace atlas::grid::global::gaussian;

namespace atlas {
namespace test {

mesh::Mesh::Ptr generate_mesh( const global::Structured& rgg )
{
  mesh::generators::Structured generate;
  generate.options.set("nb_parts",eckit::mpi::size());
  generate.options.set("part",eckit::mpi::rank());
  return mesh::Mesh::Ptr( generate( rgg ) );
}

mesh::Mesh::Ptr generate_mesh(int nlat, long lon[] )
{
  return generate_mesh( ReducedGaussian(nlat,lon) );
}


} // end namespace test
} // end namespace atlas
