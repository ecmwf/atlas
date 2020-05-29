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
#include "arrayViewHelpers.h"

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
void prff (const std::string & name, const atlas::FieldSet & s)
{
  auto & comm = mpi::comm ();
  int irank = comm.rank ();
  for (int i = 0; i < s.size (); i++)
    {
      auto v = array::make_view<T,1> (s[i]);
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

template <typename T>
void prff (const std::string & name, const atlas::Field & f)
{
  atlas::FieldSet s;
  s.add (f);
  prff<T> (name, f);
}

};

//-----------------------------------------------------------------------------

CASE( "test_gatherscatter_nflevgxngptot" ) 
{
    if (! eckit::Resource<bool> ("--nflevgxngptot", false))
      return;

    int nfield = eckit::Resource<int> ("--fields", 3);
    int nflevg = eckit::Resource<int> ("--nflevg", 10);
    atlas::StructuredGrid grid (eckit::Resource<std::string> ("--grid", "N16"));
    bool gather (eckit::Resource<bool> ("--gather", false));
    bool scatter (eckit::Resource<bool> ("--scatter", false));
    bool check (eckit::Resource<bool> ("--check", false));


    auto & comm = mpi::comm ();
    int irank = comm.rank ();
    int nproc = comm.size ();

    atlas::grid::Distribution dist (grid, atlas::grid::Partitioner ("equal_regions"));
    atlas::functionspace::StructuredColumns fs (grid, dist, atlas::util::Config ("halo", 1) 
                                              | atlas::util::Config ("periodic_points", true));
    std::vector<int> prc, ind;

    if (check)
      getprcind (fs, dist, prc, ind);

    atlas::FieldSet sloc;

    using T = long;

    const int ind_bit = 24, ind_off =  0;
    const int prc_bit = 12, prc_off = ind_off + ind_bit;
    const int fld_bit = 14, fld_off = prc_off + prc_bit;
    const int lev_bit = 14, lev_off = fld_off + fld_bit;
  
    if (check)
      {
        ATLAS_ASSERT (fs.sizeOwned () < (1 << ind_bit));
        ATLAS_ASSERT (nproc           < (1 << prc_bit));
        ATLAS_ASSERT (nfield          < (1 << fld_bit));
        ATLAS_ASSERT (nflevg          < (1 << lev_bit));
     }

    using T = long;

    auto func = [check,ind_off,prc_off,fld_off,lev_off] (int fld, int lev, int prc, int ind)
    {
      long v = 0;
      if (check)
        {
          v = (static_cast<long> (fld) << fld_off) + (static_cast<long> (lev) << lev_off) 
            + (static_cast<long> (prc) << prc_off) + (static_cast<long> (ind) << ind_off);
        }
      return v;
    };

    atlas::Field floc ("field", atlas::array::DataType::kind<T> (), {nfield, grid.size (), nflevg});

    auto v = array::make_view<T,3> (floc);

    for (int jfld = 0; jfld < nfield; jfld++)
    for (int jlev = 0; jlev < nflevg; jlev++)
    for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
      v (jfld, jloc, jlev) = func (jfld, jlev, irank, jloc);

    atlas::FieldSet sglo;

    for (int jfld = 0, count = 0; jfld < nfield; jfld++)
    for (int jlev = 0; jlev < nflevg; jlev++, count++)
      {
        int owner = count % nproc;
        std::string name = std::string ("#") + std::to_string (jfld) + std::string (",") + std::to_string (jlev);
        atlas::Field fglo = atlas::Field (name, atlas::array::DataType::kind<T> (), {irank == owner ? grid.size () : 0}); 
        fglo.metadata ().set ("owner", owner);
        sglo.add (fglo);
      }

    std::vector<ioFieldDesc> dloc;
    std::vector<ioFieldDesc> dglo;

    createIoFieldDescriptors (floc, dloc, fs.sizeOwned (), 1);
    createIoFieldDescriptors (sglo, dglo, fs.sizeOwned ());

    for (auto & d : dglo)
      d.field ().metadata ().get ("owner", d.owner ());
  
    GatherScatter gs (dist);
    gs.gather (dloc, dglo);


    for (int jfld = 0, count = 0; jfld < nfield; jfld++)
    for (int jlev = 0; jlev < nflevg; jlev++, count++)
      {
        atlas::Field fglo = sglo[count];
        int owner; fglo.metadata ().get ("owner", owner);
        if (owner == irank)
          {
            auto v = array::make_view<T,1> (fglo);
            for (int jglo = 0; jglo < grid.size (); jglo++)
              {
                T v1 = v (jglo), v2 = func (jfld, jlev, prc[jglo], ind[jglo]);
                EXPECT (v1 == v2);
              }
          }
      }

  {
  auto v = array::make_view<T,3> (floc);
  filterViewHelper<3,T>::apply (v, [](T & z){ z = 0; });

  }


}

CASE( "test_gatherscatter_ngptotxnflevg" ) 
{
    if (! eckit::Resource<bool> ("--ngptotxnflevg", false))
      return;
    int nfield = eckit::Resource<int> ("--fields", 3);
    int nflevg = eckit::Resource<int> ("--nflevg", 10);
    atlas::StructuredGrid grid (eckit::Resource<std::string> ("--grid", "N16"));
    bool gather (eckit::Resource<bool> ("--gather", false));
    bool scatter (eckit::Resource<bool> ("--scatter", false));
    bool check (eckit::Resource<bool> ("--check", false));


    auto & comm = mpi::comm ();
    int irank = comm.rank ();
    int nproc = comm.size ();

    atlas::grid::Distribution dist (grid, atlas::grid::Partitioner ("equal_regions"));
    atlas::functionspace::StructuredColumns fs (grid, dist, atlas::util::Config ("halo", 1) 
                                              | atlas::util::Config ("periodic_points", true));
    std::vector<int> prc, ind;

    if (check)
      getprcind (fs, dist, prc, ind);

    atlas::FieldSet sloc;

    using T = long;

    const int ind_bit = 24, ind_off =  0;
    const int prc_bit = 12, prc_off = ind_off + ind_bit;
    const int fld_bit = 14, fld_off = prc_off + prc_bit;
    const int lev_bit = 14, lev_off = fld_off + fld_bit;
  
    if (check)
      {
        ATLAS_ASSERT (fs.sizeOwned () < (1 << ind_bit));
        ATLAS_ASSERT (nproc           < (1 << prc_bit));
        ATLAS_ASSERT (nfield          < (1 << fld_bit));
        ATLAS_ASSERT (nflevg          < (1 << lev_bit));
     }

    using T = long;

    auto func = [check,ind_off,prc_off,fld_off,lev_off] (int fld, int lev, int prc, int ind)
    {
      long v = 0;
      if (check)
        {
          v = (static_cast<long> (fld) << fld_off) + (static_cast<long> (lev) << lev_off) 
            + (static_cast<long> (prc) << prc_off) + (static_cast<long> (ind) << ind_off);
        }
      return v;
    };

    atlas::Field floc ("field", atlas::array::DataType::kind<T> (), {nfield, nflevg, grid.size ()});


    auto v = array::make_view<T,3> (floc);

    for (int jfld = 0; jfld < nfield; jfld++)
    for (int jlev = 0; jlev < nflevg; jlev++)
    for (int jloc = 0; jloc < fs.sizeOwned (); jloc++)
      v (jfld, jlev, jloc) = func (jfld, jlev, irank, jloc);

    atlas::FieldSet sglo;

    for (int jfld = 0, count = 0; jfld < nfield; jfld++)
    for (int jlev = 0; jlev < nflevg; jlev++, count++)
      {
        int owner = count % nproc;
        std::string name = std::string ("#") + std::to_string (jfld) + std::string (",") + std::to_string (jlev);
        atlas::Field fglo = atlas::Field (name, atlas::array::DataType::kind<T> (), {irank == owner ? grid.size () : 0}); 
        fglo.metadata ().set ("owner", owner);
        sglo.add (fglo);
      }

    std::vector<ioFieldDesc> dloc;
    std::vector<ioFieldDesc> dglo;

    createIoFieldDescriptors (floc, dloc, fs.sizeOwned ());
    createIoFieldDescriptors (sglo, dglo, fs.sizeOwned ());

    for (auto & d : dglo)
      d.field ().metadata ().get ("owner", d.owner ());
  
    GatherScatter gs (dist);
    gs.gather (dloc, dglo);


    for (int jfld = 0, count = 0; jfld < nfield; jfld++)
    for (int jlev = 0; jlev < nflevg; jlev++, count++)
      {
        atlas::Field fglo = sglo[count];
        int owner; fglo.metadata ().get ("owner", owner);
        if (owner == irank)
          {
            auto v = array::make_view<T,1> (fglo);
            for (int jglo = 0; jglo < grid.size (); jglo++)
              {
                T v1 = v (jglo), v2 = func (jfld, jlev, prc[jglo], ind[jglo]);
                EXPECT (v1 == v2);
              }
          }
      }


}

CASE( "test_gatherscatter_simple" ) 
{
    if (! eckit::Resource<bool> ("--simple", false))
      return;

    int nfield = eckit::Resource<int> ("--fields", 3);
    atlas::StructuredGrid grid (eckit::Resource<std::string> ("--grid", "N16"));
    bool gather1 (eckit::Resource<bool> ("--gather1", false));
    bool gather2 (eckit::Resource<bool> ("--gather2", false));
    bool scatter1 (eckit::Resource<bool> ("--scatter1", false));
    bool scatter2 (eckit::Resource<bool> ("--scatter2", false));
    bool debug (eckit::Resource<bool> ("--debug", false));
    bool check (eckit::Resource<bool> ("--check", false));

    auto & comm = mpi::comm ();
    int irank = comm.rank ();
    int nproc = comm.size ();

    atlas::grid::Distribution dist (grid, atlas::grid::Partitioner ("equal_regions"));
    atlas::functionspace::StructuredColumns fs (grid, dist, atlas::util::Config ("halo", 1) 
                                              | atlas::util::Config ("periodic_points", true));

    std::vector<int> prc, ind;

    if (check)
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
        ATLAS_ASSERT (nfield          < (1 << fld_bit));
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
      ATLAS_TRACE_SCOPE ("checkgather")
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

    auto checkscatter = [fs, func, irank] (const atlas::FieldSet & sloc)
    {
      ATLAS_TRACE_SCOPE ("checkscatter")
      {
      for (int i = 0; i < sloc.size (); i++)
        {
          const auto & floc = sloc[i];
          const auto v = array::make_view<T,1> (floc);
          for (int j = 0; j < fs.sizeOwned (); j++)
            {
              T v1 = v[j], v2 = func (i, irank, j);
              EXPECT_EQ (v1, v2);
            }
        }
      }
    };

    for (int i = 0; i < nfield; i++)
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

        if (scatter1)
          {
            ATLAS_TRACE_SCOPE ("test_scatter1")
            {
            for (int i = 0; i < sloc.size (); i++)
              {
                auto v = array::make_view<T,1> (sloc[i]);
                for (int j = 0; j < fs.sizeOwned (); j++)
                  v (j) = 0;
              }

            fs.scatter (sglo1, sloc);

            if (check) checkscatter (sloc);
            }
          }
      }

    if (gather2)
      {
        GatherScatter gs (dist);
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

          if (scatter2)
            {
              ATLAS_TRACE_SCOPE ("test_scatter2")
              {
              for (int i = 0; i < sloc.size (); i++)
                {
                  auto v = array::make_view<T,1> (sloc[i]);
                  for (int j = 0; j < fs.sizeOwned (); j++)
                    v (j) = 0;
                }

              gs.scatter (dglo, dloc);

              if (check) checkscatter (sloc);
              }
            }
        }
      }


  // Prevent error in Atlas barrier
  comm.barrier ();

}


}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
