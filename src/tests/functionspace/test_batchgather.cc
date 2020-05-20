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
#include "atlas/runtime/Trace.h"

#include "GatherScatter.h"
#include "ioFieldDesc.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

namespace
{

void getprcind (const atlas::functionspace::StructuredColumns & fs, const atlas::grid::Distribution & dist,
                std::vector<int> & prc, std::vector<int> & ind)
{
  auto & comm = mpi::comm ();
  auto & grid = fs.grid ();
  prc.resize (grid.size ());
  ind.resize (grid.size ());

  for (int i = 0; i < grid.size (); i++)
    prc[i] = dist.partition (i);
 
  {
    atlas::Field indloc = fs.createField<int> (atlas::util::Config ("name", "ind") 
                                             | atlas::util::Config ("owner", 0));
    atlas::Field indglo ("ind", &ind[0], {grid.size ()});
    auto v = array::make_view<int,1> (indloc);
    for (int i = 0; i < fs.sizeOwned (); i++)
      v (i) = i;
    fs.gather (indloc, indglo);
    comm.broadcast (ind, 0);
  }
}

template <typename T>
void prff (const std::string & name, const atlas::FieldSet & sglo)
{
  auto & comm = mpi::comm ();
  int irank = comm.rank ();
  for (int i = 0; i < sglo.size (); i++)
    {
      auto v = array::make_view<T,1> (sglo[i]);
      if (v.size () == 0)
        continue;
  
      char tmp[128];
      sprintf (tmp, "%s.%8.8d.%8.8d.txt", name.c_str (), irank, i);
      FILE * fp = fopen (tmp, "w");
      for (int i = 0; i < v.size (); i++)
        fprintf (fp, "%8d > %8d\n", i, v (i));
      fclose (fp);
    }
}

};

//-----------------------------------------------------------------------------

CASE( "test_gatherscatter" ) 
{
    int nfields = eckit::Resource<int> ("--fields", 3);
    atlas::StructuredGrid grid (eckit::Resource<std::string> ("--grid", "N16"));
    bool gather1 (eckit::Resource<bool> ("--gather1", false));
    bool gather2 (eckit::Resource<bool> ("--gather2", false));
    bool debug (eckit::Resource<bool> ("--debug", false));
    bool check (eckit::Resource<bool> ("--check", false));

    auto & comm = mpi::comm ();
    int irank = comm.rank ();
    int nproc = comm.size ();

    atlas::grid::Distribution dist (grid, atlas::grid::Partitioner ("equal_regions"));
    atlas::functionspace::StructuredColumns fs (grid, dist, atlas::util::Config ("halo", 1) 
                                              | atlas::util::Config ("periodic_points", true));

    std::vector<int> prc, ind;

    getprcind (fs, dist, prc, ind);

    atlas::FieldSet sloc;
    atlas::FieldSet sglo1;
    atlas::FieldSet sglo2;

    using T = long;

    const int ind_bit = 21, ind_off =  0;
    const int prc_bit = 21, prc_off = ind_off + ind_bit;
    const int fld_bit = 21, fld_off = prc_off + prc_bit;
  
    if (check)
      {
        ATLAS_ASSERT (fs.sizeOwned () < (1 << ind_bit));
        ATLAS_ASSERT (nproc           < (1 << prc_bit));
        ATLAS_ASSERT (nfields         < (1 << fld_bit));
     }

    auto func = [check,ind_off,prc_off,fld_off] (int fld, int prc, int ind)
    {
      long v = 0;
      if (check)
        {
          v = (static_cast<long> (fld) << fld_off) 
            + (static_cast<long> (prc) << prc_off) 
            + (static_cast<long> (ind) << ind_off);
        }
      return v;
    };

    auto checkgather = [ind, prc, grid, func, irank] (const atlas::FieldSet & sglo)
    {
      ATLAS_TRACE_SCOPE ("check")
      {
      for (int i = 0; i < sglo.size (); i++)
        {
          const auto & fglo = sglo[i];
          int owner;
          EXPECT (fglo.metadata ().get ("owner", owner));
          if (owner == irank)
            {
              const auto v = array::make_view<T,1> (fglo);
              for (int j = 0; j < grid.size (); j++)
                {
                  T v1 = v[j], v2 = func (i, prc[j], ind[j]);
                  EXPECT_EQ (v1, v2);
                }
            }
        }
      }
    };

    for (int i = 0; i < nfields; i++)
      {
        int owner = i % nproc;

        std::string name = std::string ("#") + std::to_string (i);
        atlas::Field floc = fs.createField<T> (atlas::util::Config ("name", name) 
                                            |  atlas::util::Config ("owner", owner));
        auto v = array::make_view<T,1> (floc);
        for (int j = 0; j < fs.sizeOwned (); j++)
          v[j] = func (i, irank, j);
        sloc.add (floc);

        if (gather1)
          {
            Field fglo1 = Field (name, atlas::array::DataType::kind<T> (), {irank == owner ? grid.size () : 0}); 
            fglo1.metadata ().set ("owner", owner);
            sglo1.add (fglo1);
          }

        if (gather2)
          {
            Field fglo2 = Field (name, atlas::array::DataType::kind<T> (), {irank == owner ? grid.size () : 0}); 
            fglo2.metadata ().set ("owner", owner);
            sglo2.add (fglo2);
          }

      }

    if (gather1)
      {
        {
          atlas::FieldSet loc, glo;
          fs.gather (loc, glo);
        }
        ATLAS_TRACE_SCOPE ("test_gather1")
        {
          fs.gather (sloc, sglo1);
          if (check) checkgather (sglo1);
          if (debug) prff<T> ("sglo1", sglo1);
        }
      }

    if (gather2)
      {
        GatherScatter gs (grid, dist);
        ATLAS_TRACE_SCOPE ("test_gather2")
        {
          std::vector<ioFieldDesc> dloc;
          std::vector<ioFieldDesc> dglo;

          ATLAS_TRACE_SCOPE ("create io descriptors")
          {
          createIoFieldDescriptors (sloc,  dloc, fs.sizeOwned ());
          createIoFieldDescriptors (sglo2, dglo, grid.size ());
          for (auto & d : dglo)
            d.field ().metadata ().get ("owner", d.owner ());
          }

          gs.gather (dloc, dglo);

          if (debug) prff<T> ("sglo2", sglo2);
          if (check) checkgather (sglo2);
        }
      }

  // Avoid Atlas barrier
  comm.barrier ();

}


}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
