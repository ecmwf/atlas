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


#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>



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

CASE( "test_broken" ) 
{
  auto & comm = atlas::mpi::comm ();
  int nproc = comm.size ();

  const int nx = 400, ny = 200;

  StructuredGrid grid = Grid (std::string ("L") + std::to_string (nx) + "x" + std::to_string (ny));

  const int icnt[] = {26800, 26800, 26400};
  int partition[grid.size ()];

  for (int irank = 0, ioff = 0; irank < 3; irank++)
    {
      for (int jglo = ioff; jglo < ioff + icnt[irank]; jglo++)
        partition[jglo] = irank;
      ioff += icnt[irank];
    }

  printf (" icnt = %8d\n", icnt[0]+icnt[1]+icnt[2]);
  printf (" grid.size () = %8d\n", grid.size ());

  EXPECT (grid.size () == icnt[0]+icnt[1]+icnt[2]);
  EXPECT (comm.size () == 3);
  
  atlas::grid::Distribution dist (comm.size (), grid.size (), partition);
  atlas::functionspace::StructuredColumns fs (grid, dist, atlas::util::Config ("halo", 1) | Config ("periodic_points", true));

  auto ij = getIJ (fs);
  fs.haloExchange (ij);

}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
