/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/parallel/mpi/mpi.h"
#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_functionspace_StructuredColumns_batchgather" ) {
    int nfields = eckit::Resource<int> ("--fields", 3);
    atlas::StructuredGrid grid (eckit::Resource<std::string> ("--grid", "N16"));

    auto & comm = mpi::comm ();
    int irank = comm.rank ();
    int nproc = comm.size ();

    atlas::grid::Distribution dist (grid, atlas::grid::Partitioner ("equal_regions"));
    atlas::functionspace::StructuredColumns fs (grid, dist, atlas::util::Config ("halo", 1) 
                                              | atlas::util::Config ("periodic_points", true));

    std::vector<int> prc (grid.size ());
    std::vector<int> ind (grid.size ());

    for (int i = 0; i < grid.size (); i++)
      prc[i] = dist.partition (i);
   
    {
      atlas::Field indloc = fs.createField<int> (atlas::util::Config ("name", "ind") | atlas::util::Config ("owner", 0));
      atlas::Field indglo ("ind", &ind[0], {grid.size ()});
      auto v = array::make_view<int,1> (indloc);
      for (int i = 0; i < fs.sizeOwned (); i++)
        v (i) = i;
      fs.gather (indloc, indglo);
      comm.broadcast (ind, 0);
    }


    atlas::FieldSet sloc;
    atlas::FieldSet sglo;

    auto func = [] (int fld, int prc, int ind)
    {
      double v = 1000.0 * static_cast<double> (fld) + prc + static_cast<double> (ind) / 1000.0;
      return v;
    };

    for (int i = 0; i < nfields; i++)
      {
        int owner = i % nproc;
        std::string name = std::string ("#") + std::to_string (i);
        atlas::Field floc = fs.createField<double> (atlas::util::Config ("name", name) 
                                                  | atlas::util::Config ("owner", owner));
        auto v = array::make_view<double,1> (floc);
        for (int j = 0; j < fs.sizeOwned (); j++)
          v[j] = func (i, irank, j);
        sloc.add (floc);
        Field fglo = Field (name, atlas::array::DataType::kind<double> (), {grid.size ()}); 
        fglo.metadata ().set ("owner", owner);
        sglo.add (fglo);
      }

    fs.gather (sloc, sglo);

    for (int i = 0; i < nfields; i++)
      {
        const auto & fglo = sglo[i];
        int owner;
        EXPECT (fglo.metadata ().get ("owner", owner));
        if (owner == irank)
          {
            const auto v = array::make_view<double,1> (fglo);
            for (int j = 0; j < grid.size (); j++)
              {
                double v1 = v[j], v2 = func (i, prc[j], ind[j]);
                EXPECT_EQ (v1, v2);
              }
          }
      }


}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
