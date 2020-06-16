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

CASE( "test_gatherscatter_NPROMAxNFLEVGxNGPBLKS" ) 
{
    if (! eckit::Resource<bool> ("--NPROMAxNFLEVGxNGPBLKS", false))
      return;

    int nproma = eckit::Resource<int> ("--nproma", 24);
    int nflevg = eckit::Resource<int> ("--nflevg", 10);
    bool gather (eckit::Resource<bool> ("--gather", false));
    bool scatter (eckit::Resource<bool> ("--scatter", false));
    bool check (eckit::Resource<bool> ("--check", false));

    atlas::StructuredGrid grid (eckit::Resource<std::string> ("--grid", "N16"));

    auto & comm = mpi::comm ();
    int irank = comm.rank ();
    int nproc = comm.size ();

    atlas::grid::Distribution dist (grid, atlas::grid::Partitioner ("equal_regions"));
    atlas::functionspace::StructuredColumns fs (grid, dist, atlas::util::Config ("halo", 1) 
                                              | atlas::util::Config ("periodic_points", true));

    int ngptot = fs.sizeOwned ();
    int ngpblks = ngptot / nproma;

    if (fs.sizeOwned () % nproma)
      ngpblks++;

    std::vector<int> prc, ind;

    if (check)
      getprcind (fs, dist, prc, ind);

    atlas::FieldSet sloc;

    using T = long;

    const int ind_bit = 21, ind_off =  0;
    const int prc_bit = 21, prc_off = ind_off + ind_bit;
    const int lev_bit = 22, lev_off = prc_off + prc_bit;
  
    if (check)
      {
        ATLAS_ASSERT (fs.sizeOwned () < (1 << ind_bit));
        ATLAS_ASSERT (nproc           < (1 << prc_bit));
        ATLAS_ASSERT (nflevg          < (1 << lev_bit));
     }

    using T = long;

    auto func = [check,ind_off,prc_off,lev_off] (int lev, int prc, int ind)
    {
      long v = 0;
      v = (static_cast<long> (lev) << lev_off) 
        + (static_cast<long> (prc) << prc_off) 
        + (static_cast<long> (ind) << ind_off);
      return v;
    };

    // Distributed field, Fortran dimensions (1:NPROMA,1:NFLEVG,1:NGPBLKS)
    atlas::Field floc ("field", atlas::array::DataType::kind<T> (), {ngpblks, nflevg, nproma});

    // Gather our multi-field, multi-level Atlas field to a set of fields (1:NGPTOTG)
    atlas::FieldSet sglo;

    for (int jlev = 0; jlev < nflevg; jlev++)
      {
        int owner = jlev % nproc;  // RR distribution
        std::string name = std::string ("#") + std::to_string (jlev);
        atlas::Field fglo = atlas::Field (name, atlas::array::DataType::kind<T> (), {irank == owner ? grid.size () : 0}); 
        fglo.metadata ().set ("owner", owner);
        sglo.add (fglo);
      }

    // IO descriptors
    std::vector<ioFieldDesc> dloc;
    std::vector<ioFieldDesc> dglo;

    createIoFieldDescriptors (sglo, dglo); // Default for grid dimension is inner dimension
    createIoFieldDescriptorsBlocked (floc, dloc, 0, 2, fs.sizeOwned ()); // NGPBLKS dimension, NPROMA dimension, NGPTOT

    // Set target processor for all fields
    for (auto & d : dglo)
      d.field ().metadata ().get ("owner", d.owner ());

    GatherScatter gs (dist);


    auto walkLoc = [&floc, ngpblks, nproma, nflevg, ngptot, func, irank] (std::function<void(T &, T)> f)
    {
      auto vloc = array::make_view<T,3> (floc);
      viewLoop (vloc, 
                [func,f,&irank,&nproma] (T & x, int jblk, int jlev, int jlon) 
                { f (x, func (jlev, irank, jlon + jblk * nproma)); });
std::cout << " walkLoc " << std::endl;
    };

    auto walkGlo = [&sglo, nflevg, &grid, &irank, &prc, &ind, func] (std::function<void(T &, T)> f)
    {
      for (int jlev = 0; jlev < nflevg; jlev++)
        {
          auto fglo = sglo[jlev];
          int owner;
          fglo.metadata ().get ("owner", owner);
          if (irank == owner)
            {
              auto vglo = array::make_view<T,1> (fglo);
              for (int jglo = 0; jglo < grid.size (); jglo++)
                f (vglo (jglo), func (jlev, prc[jglo], ind[jglo]));
            }
        }
    };

    auto set   = [](T & v, T r) { v = r; };
    auto cmp   = [](T & v, T r) { EXPECT (v == r); };
    auto clear = [](T & v, T r) { v = 0; };

    if (gather)
      {
        if (check)
          {
            walkLoc (set);
            walkGlo (clear);
          }
        gs.gather (dloc, dglo);
        if (check)
          walkGlo (cmp);
      }

   if (scatter)
     {
       if (check)
         {
           walkGlo (set);
           walkLoc (clear);
         }
       gs.scatter (dglo, dloc);
       if (check)
         walkGlo (cmp);
     }
}

CASE( "test_gatherscatter_NFLEVGxNGPTOT" ) 
{
    if (! eckit::Resource<bool> ("--NFLEVGxNGPTOT", false))
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
      v = (static_cast<long> (fld) << fld_off) + (static_cast<long> (lev) << lev_off) 
        + (static_cast<long> (prc) << prc_off) + (static_cast<long> (ind) << ind_off);
      return v;
    };


    // Distributed field, Fortran dimensions (1:NFLEVG,1:NGPTOT,1:NFIELDS)
    atlas::Field floc ("field", atlas::array::DataType::kind<T> (), {nfield, fs.size (), nflevg});
    atlas::FieldSet sglo;

    for (int jfld = 0, count = 0; jfld < nfield; jfld++)
    for (int jlev = 0; jlev < nflevg; jlev++, count++)
      {
        int owner = count % nproc;  // RR distribution
        std::string name = std::string ("#") + std::to_string (jfld) + std::string (",") + std::to_string (jlev);
        atlas::Field fglo = atlas::Field (name, atlas::array::DataType::kind<T> (), {irank == owner ? grid.size () : 0}); 
        fglo.metadata ().set ("owner", owner);
        sglo.add (fglo);
      }

    auto walkLoc = [&floc, nfield, nflevg, &fs, irank, func] (std::function<void(T &, T)> f)
    {
      auto vloc = array::make_view<T,3> (floc);
      viewLoop (vloc, [func,&irank,f] (T & x, int jfld, int jloc, int jlev) { f (x, func (jfld, jlev, irank, jloc)); });
    };

    auto walkGlo = [&sglo, nfield, nflevg, irank, func, &grid, &prc, &ind] (std::function<void(T &, T)> f)
    {
      for (int jfld = 0, count = 0; jfld < nfield; jfld++)
      for (int jlev = 0; jlev < nflevg; jlev++, count++)
        {
          atlas::Field fglo = sglo[count];
          int owner; EXPECT (fglo.metadata ().get ("owner", owner));
          if (owner == irank)
            {
              auto v = array::make_view<T,1> (fglo);
              for (int jglo = 0; jglo < grid.size (); jglo++)
                f (v (jglo), func (jfld, jlev, prc[jglo], ind[jglo]));
            }
        }
    };

    auto set   = [](T & v, T r) { v = r; };
    auto cmp   = [](T & v, T r) { EXPECT (v == r); };
    auto clear = [](T & v, T r) { v = 0; };

    // IO descriptors
    std::vector<ioFieldDesc> dloc;
    std::vector<ioFieldDesc> dglo;

    createIoFieldDescriptors (floc, dloc, fs.sizeOwned (), 1); // Grid dimension is 1
    createIoFieldDescriptors (sglo, dglo);    // Default for grid dimension is inner dimension

    // Set target processor for all fields
    for (auto & d : dglo)
      d.field ().metadata ().get ("owner", d.owner ());
  
    GatherScatter gs (dist);

    // Gather

    if (gather)
      {
        if (check)
          {
            walkLoc (set); 
            walkGlo (clear);
          }

        gs.gather (dloc, dglo);

        if (check)
          walkGlo (cmp);
      }

  // Scatter

  if (scatter)
    {
      if (check)
        {
          walkGlo (set);
          walkLoc (clear);
        }

      gs.scatter (dglo, dloc);

      if (check)
        walkLoc (cmp);
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
    atlas::FieldSet sglo;

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
      v = (static_cast<long> (fld) << fld_off) 
        + (static_cast<long> (prc) << prc_off) 
        + (static_cast<long> (ind) << ind_off);
      return v;
    };

    auto walkGlo = [&sglo, &ind, &prc, &grid, func, irank] (std::function<void(T &, T)> f)
    {
      for (int i = 0; i < sglo.size (); i++)
        {
          auto & fglo = sglo[i];
          int owner;
          EXPECT (fglo.metadata ().get ("owner", owner));
          if (owner == irank)
            {
              auto v = array::make_view<T,1> (fglo);
              for (int j = 0; j < grid.size (); j++)
                f (v[j], func (i, prc[j], ind[j]));
            }
        }
    };

    auto walkLoc = [&sloc, &fs, func, irank] (std::function<void(T &, T)> f)
    {
      for (int i = 0; i < sloc.size (); i++)
        {
          auto & floc = sloc[i];
          auto v = array::make_view<T,1> (floc);
          for (int j = 0; j < fs.sizeOwned (); j++)
            f (v[j], func (i, irank, j));
        }
    };


    for (int i = 0; i < nfield; i++)
      {
        int owner = i % nproc;

        std::string name = std::string ("#") + std::to_string (i);
        atlas::Field floc = fs.createField<T> (atlas::util::Config ("name", name) 
                                            |  atlas::util::Config ("owner", owner));

        sloc.add (floc);

        Field fglo = Field (name, atlas::array::DataType::kind<T> (), {irank == owner ? grid.size () : 0}); 
        fglo.metadata ().set ("owner", owner);
        sglo.add (fglo);
      }

    auto set   = [](T & v, T r) { v = r; };
    auto cmp   = [](T & v, T r) { EXPECT (v == r); };
    auto clear = [](T & v, T r) { v = 0; };

    if (gather1)
      {
        if (check)
          {
            walkLoc (set);
            walkGlo (clear);
          }

        ATLAS_TRACE_SCOPE ("test_gather1")
        {
          fs.gather (sloc, sglo);
        }

        if (check) 
          walkGlo (cmp);

        if (debug) prff<T> ("sglo1", sglo);
     }

    if (scatter1)
      {
        if (check)
          {
            walkGlo (set);
            walkLoc (clear);
          }
        ATLAS_TRACE_SCOPE ("test_scatter1")
        {
          fs.scatter (sglo, sloc);
        };
        if (check) 
          walkLoc (cmp);
      }

    if (gather2)
      {
        GatherScatter gs (dist);
        if (check)
          {
            walkLoc (set);
            walkGlo (clear);
          }
        ATLAS_TRACE_SCOPE ("test_gather2")
        {
          std::vector<ioFieldDesc> dloc;
          std::vector<ioFieldDesc> dglo;

          ATLAS_TRACE_SCOPE ("create io descriptors")
          {
            createIoFieldDescriptors (sloc,  dloc, fs.sizeOwned ());
            createIoFieldDescriptors (sglo, dglo);
            for (auto & d : dglo)
              EXPECT (d.field ().metadata ().get ("owner", d.owner ()));
          }

          gs.gather (dloc, dglo);
        };

        if (check) 
          walkGlo (cmp);

        if (debug) prff<T> ("sglo2", sglo);
      }

    if (scatter2)
      {
        GatherScatter gs (dist);
        if (check)
          {
            walkGlo (set);
            walkLoc (clear);
          }
        ATLAS_TRACE_SCOPE ("test_scatter2")
        {
          std::vector<ioFieldDesc> dloc;
          std::vector<ioFieldDesc> dglo;

          ATLAS_TRACE_SCOPE ("create io descriptors")
          {
            createIoFieldDescriptors (sloc,  dloc, fs.sizeOwned ());
            createIoFieldDescriptors (sglo, dglo);
            for (auto & d : dglo)
              EXPECT (d.field ().metadata ().get ("owner", d.owner ()));
          }

          gs.scatter (dglo, dloc);
        };

        if (check) 
          walkLoc (cmp);
      }


  // Prevent error in Atlas barrier
  comm.barrier ();

}


}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
