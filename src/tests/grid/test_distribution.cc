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
#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"


#include "tests/AtlasTestEnvironment.h"

using Grid   = atlas::Grid;
using Config = atlas::util::Config;

namespace atlas {
namespace test {

atlas::FieldSet 
getIJ (const atlas::functionspace::StructuredColumns & fs) 
{
  atlas::FieldSet ij;
  
  auto i = atlas::array::make_view<int,1> (fs.index_i ());
  auto j = atlas::array::make_view<int,1> (fs.index_j ());

  ij.add (fs.index_i ());
  ij.add (fs.index_j ());

  return ij;
}


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

  bool same = (ny % nproc) == 0; // Compare light & regular distributions when possible

  if (same)
  for (int i = 0; i < grid.size (); i++)
    EXPECT (dist1.partition (i) == dist2.partition (i));

  for (int iy = 0, jglo = 0; iy < ny; iy++)
    {
      int jglo0 = jglo;
      for (int ix = 0; ix < nx; ix++, jglo++)
        EXPECT (dist2.partition (jglo) == dist2.partition (jglo0));
    }


  atlas::functionspace::StructuredColumns fs1 (grid, dist1, atlas::util::Config ("halo", 1) | Config ("periodic_points", true));
  atlas::functionspace::StructuredColumns fs2 (grid, dist2, atlas::util::Config ("halo", 1) | Config ("periodic_points", true));

  auto ij1 = getIJ (fs1);
  auto ij2 = getIJ (fs2);

  fs1.haloExchange (ij1);
  fs2.haloExchange (ij2);

  if (same)
    {
      EXPECT (fs1.size () == fs1.size ());
      EXPECT (fs1.sizeOwned () == fs1.sizeOwned ());

      auto i1 = atlas::array::make_view<int,1> (ij1[0]);
      auto j1 = atlas::array::make_view<int,1> (ij1[1]);
      auto i2 = atlas::array::make_view<int,1> (ij2[0]);
      auto j2 = atlas::array::make_view<int,1> (ij2[1]);

      for (int k = 0; k < fs1.sizeOwned (); k++)
        {
          EXPECT (i1[k] == i2[k]);
          EXPECT (j1[k] == j2[k]);
        }
    }


  for (int j = fs2.j_begin_halo (); j < fs2.j_end_halo (); j++)
     {
       EXPECT (fs2.i_begin_halo (j) == -1);
       EXPECT (fs2.i_end_halo (j) == nx + 2);
     }

}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
