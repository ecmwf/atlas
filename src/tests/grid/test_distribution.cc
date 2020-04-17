/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <iomanip>
#include <sstream>

#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_nbands" ) 
{
  auto & comm = atlas::mpi::comm ();
  int nproc = comm.size ();

  StructuredGrid grid = Grid ("L200x100");
  atlas::grid::Distribution dist (grid, atlas::util::Config ("type", "checkerboard") | Config ("nbands", nproc));

  for (int i = 1; i < grid.size (); i++)
    EXPECT (dist.partition (i-1) <= dist.partition (i));
}

CASE( "test_light" ) 
{
  auto & comm = atlas::mpi::comm ();
  int nproc = comm.size ();

  const int nx = 400, ny = 200;

  StructuredGrid grid = Grid (std::string ("L") + std::to_string (nx) + "x" + std::to_string (ny));
  atlas::grid::Distribution dist1 (grid, atlas::util::Config ("type", "checkerboard") | Config ("nbands", nproc));

  atlas::grid::Distribution dist2 (grid, Config ("light", true) | Config ("blocksize", nx));
  
  EXPECT (dist2.footprint () < 100);

  if ((ny % nproc) == 0)
  for (int i = 0; i < grid.size (); i++)
    EXPECT (dist1.partition (i) == dist2.partition (i));

  for (int iy = 0, jglo = 0; iy < ny; iy++)
    {
      int jglo0 = jglo;
      for (int ix = 0; ix < nx; ix++, jglo++)
        EXPECT (dist2.partition (jglo) == dist2.partition (jglo0));
    }

}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
