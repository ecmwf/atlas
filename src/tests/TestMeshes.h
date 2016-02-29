/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/atlas_config.h"
#include "atlas/grid/GaussianLatitudes.h"
#include "atlas/grid/grids.h"
#include "atlas/mesh/generators/ReducedGridMeshGenerator.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/util/parallel/mpi/mpi.h"

using namespace atlas;
using namespace atlas::grid;

namespace atlas {
namespace test {

class TestGrid: public ReducedGrid {
public:
  TestGrid(int N, long lon[]);
};

TestGrid::TestGrid(int N, long lon[])
{
  std::vector<double> lats(N);
  grid::gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon,internals::DEG);
}

mesh::Mesh::Ptr generate_mesh( const ReducedGrid& rgg )
{
  mesh::generators::ReducedGridMeshGenerator generate;
  generate.options.set("nb_parts",eckit::mpi::size());
  generate.options.set("part",eckit::mpi::rank());
  return mesh::Mesh::Ptr( generate( rgg ) );
}

mesh::Mesh::Ptr generate_mesh(int nlat, long lon[] )
{
  return generate_mesh( TestGrid(nlat,lon) );
}


} // end namespace test
} // end namespace atlas
