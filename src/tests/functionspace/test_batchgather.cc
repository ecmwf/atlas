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

template <typename T, typename I>
T reorder (const T & vec, const I & ord)
{
  T v;
  v.reserve (ord.size ());
  for (typename I::value_type i = 0; i < ord.size (); i++)
    v.push_back (vec[ord[i]]);
  return v;
}

template <typename T>
void integrate (T & v)
{
  v[0].off = 0;
  for (size_t i = 1; i < v.size (); i++)
    v[i].off = v[i-1].off + v[i-1].len;
}

class mapping_t
{
private:

  class locprc_t
  {
  public:
    int loc = std::numeric_limits<int>::min ();
    int prc = std::numeric_limits<int>::min ();
  };

  int max, nprc;

  std::vector<int> _prcloc2glo;
  std::vector<locprc_t> _glo2prcloc;

  const atlas::StructuredGrid & grid;
  const atlas::grid::Distribution & dist;

public:
  mapping_t (const atlas::StructuredGrid & _grid, const atlas::grid::Distribution & _dist);

  int prcloc2glo (int iprc, int jloc) const
  {
    return _prcloc2glo[iprc * max + jloc];
  }

  const locprc_t & glo2prcloc (int jglo) const
  {
    return _glo2prcloc[jglo];
  }
 

};

mapping_t::mapping_t (const atlas::StructuredGrid & _grid, const atlas::grid::Distribution & _dist)
 : grid (_grid), dist (_dist)
{
  nprc = dist.nb_partitions ();

  std::vector<int> count (nprc);
  std::vector<int> ind (grid.size ());

  for (int i = 0; i < grid.size (); i++)
    ind[i] = count[dist.partition (i)]++;

  max = *std::max_element (count.begin (), count.end ());

  _prcloc2glo.resize (max * nprc, std::numeric_limits<int>::min ());
  _glo2prcloc.resize (grid.size ());

  for (int i = 0; i < grid.size (); i++)
    {
      _prcloc2glo[max * dist.partition (i) + ind[i]] = i;
      _glo2prcloc[i].loc = ind[i];
      _glo2prcloc[i].prc = dist.partition (i);
    }

}

void gather (const atlas::StructuredGrid & grid, const atlas::grid::Distribution & dist, 
             std::vector<ioFieldDesc> & floc, std::vector<ioFieldDesc> & fglo)
{
  ATLAS_ASSERT (floc.size () == fglo.size ());
  size_t nfld = floc.size ();

  auto & comm = mpi::comm ();
  int nprc = comm.size ();
  int lprc = comm.rank ();

  class offlen_t
  {
  public:
    size_t off = 0, len = 0;
  };

  mapping_t mapping (grid, dist);

  size_t ldim = dist.nb_pts ()[lprc];

  // Sort fields by owner
  
  std::vector<int> isort (nfld);

  std::iota (std::begin (isort), std::end (isort), 0);
  std::sort (std::begin (isort), std::end (isort), 
             [&fglo] (int a, int b) { return fglo[a].owner () < fglo[b].owner (); });

  floc  = reorder (floc,  isort);
  fglo  = reorder (fglo,  isort);

  // SEND

  std::vector<offlen_t> fld_send (nfld + 1);
  std::vector<offlen_t> prc_send (nprc + 1);

  for (int jfld = 0; jfld < nfld; jfld++)
    {
      int owner = fglo[jfld].owner ();
      fld_send[jfld].len = floc[jfld].size ();
      prc_send[owner].len += floc[jfld].size ();
    }

printf (" dist.nb_pts ()[lprc] = %8d\n", dist.nb_pts ()[lprc]);

  integrate (fld_send);
  integrate (prc_send);

  // Pack send buffer

  std::vector<byte> buf_send (fld_send.back ().off);

  for (int jfld = 0; jfld < nfld; jfld++)
    floc[jfld].pack (&buf_send[fld_send[jfld].off]);

  // RECV

  std::vector<offlen_t> prc_recv (nprc + 1);
  std::vector<offlen_t> fld_recv (nfld + 1);

  for (int jfld = 0; jfld < nfld; jfld++)
    if (lprc == fglo[jfld].owner ())
      fld_recv[jfld].len = fglo[jfld].dlen ();

  integrate (fld_recv);

  for (int iprc = 0; iprc < nprc; iprc++)
     prc_recv[iprc].len = dist.nb_pts ()[iprc] * fld_recv.back ().off;

  integrate (prc_recv);

  std::vector<byte> buf_recv (prc_recv.back ().off);

  std::vector<eckit::mpi::Request> rqr;

  for (int iprc = 0; iprc < nprc; iprc++)
    if (prc_recv[iprc].len > 0)
      rqr.push_back (comm.iReceive (&buf_recv[prc_recv[iprc].off], 
                                    prc_recv[iprc].len, iprc, 100));

  comm.barrier ();

  std::vector<eckit::mpi::Request> rqs;

  for (int iprc = 0; iprc < nprc; iprc++)
    if (prc_send[iprc].len > 0)
      rqs.push_back (comm.iSend (&buf_send[prc_send[iprc].off], 
                                 prc_send[iprc].len, iprc, 100));

  for (auto & r : rqr)
    comm.wait (r);

  for (auto & r : rqs)
    comm.wait (r);

  // Unpack RECV buffer

  for (int iprc = 0; iprc < nprc; iprc++)
    if (prc_recv[iprc].len > 0)
      {
        int off = prc_recv[iprc].off;
        for (int jfld = 0; jfld < nfld; jfld++)
          {
            if (fld_recv[jfld].len > 0)
              {
                for (int jloc = 0, k = 0; jloc < dist.nb_pts ()[iprc]; jloc++)
                for (int j = 0; j < fld_recv[jfld].len; j++, k++)
                  fglo[jfld](mapping.prcloc2glo (iprc, jloc), j) = buf_recv[off+k];
                off += dist.nb_pts ()[iprc] * fglo[jfld].dlen ();
              }
          }
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

CASE( "test_functionspace_StructuredColumns_batchgather" ) {
    int nfields = eckit::Resource<int> ("--fields", 3);
    atlas::StructuredGrid grid (eckit::Resource<std::string> ("--grid", "N16"));

    auto & comm = mpi::comm ();
    int irank = comm.rank ();
    int nproc = comm.size ();

    atlas::grid::Distribution dist (grid, atlas::grid::Partitioner ("equal_regions"));
    atlas::functionspace::StructuredColumns fs (grid, dist, atlas::util::Config ("halo", 1) 
                                              | atlas::util::Config ("periodic_points", true));

    std::vector<int> prc, ind;

    getprcind (fs, dist, prc, ind);

    atlas::FieldSet sloc;
    atlas::FieldSet sglo1, sglo2;

    {
    char tmp[128];
    sprintf (tmp, "prcind.%8.8d.txt", irank);
    FILE * fp = fopen (tmp, "w");
    for (int i = 0; i < grid.size (); i++)
      fprintf (fp, " %8d > %8d, %8d\n", i, prc[i], ind[i]);
    fclose (fp);
    }

    using T = long;

    const int ind_bit = 21, ind_off =  0;
    const int prc_bit = 21, prc_off = ind_off + ind_bit;
    const int fld_bit = 21, fld_off = prc_off + prc_bit;
  
    ATLAS_ASSERT (fs.sizeOwned () < (1 << ind_bit));
    ATLAS_ASSERT (nproc           < (1 << prc_bit));
    ATLAS_ASSERT (nfields         < (1 << fld_bit));

    auto func = [&ind_off,&prc_off,&fld_off] (int fld, int prc, int ind)
    {
      long v = (static_cast<long> (fld) << fld_off) 
             + (static_cast<long> (prc) << prc_off) 
             + (static_cast<long> (ind) << ind_off);
      return v;
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

        Field fglo = Field (name, atlas::array::DataType::kind<T> (), {irank == owner ? grid.size () : 0}); 
        fglo.metadata ().set ("owner", owner);
        sglo1.add (fglo);

        Field fglo2 = Field (name, atlas::array::DataType::kind<T> (), {irank == owner ? grid.size () : 0}); 
        fglo2.metadata ().set ("owner", owner);
        sglo2.add (fglo2);

      }


    fs.gather (sloc, sglo1);

    {
      std::vector<ioFieldDesc> dloc;
      std::vector<ioFieldDesc> dglo;

      createIoFieldDescriptors (sloc,  dloc, fs.sizeOwned ());
      createIoFieldDescriptors (sglo2, dglo, grid.size ());
     
      for (auto & d : dglo)
        d.field ().metadata ().get ("owner", d.owner ());

      printf (" dloc.size () = %8d\n", dloc.size ());
      printf (" dglo.size () = %8d\n", dglo.size ());

      gather (grid, dist, dloc, dglo);


      prff<T> ("sglo2", sglo2);

    }


    prff<T> ("sglo1", sglo1);


    auto cmp = [ind, prc, grid, func, irank] (const atlas::FieldSet & sglo)
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
    };

    cmp (sglo1);
    cmp (sglo2);

}


}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
