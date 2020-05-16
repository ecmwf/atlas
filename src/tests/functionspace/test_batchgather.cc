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

namespace
{


typedef char byte;

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
  T v (ord.size ());
  for (typename I::value_type i = 0; i < ord.size (); i++)
    v[i] = vec[ord[i]];
  return v;
}

template <typename T>
void integrate (T & v)
{
  v[0].off = 0;
  for (size_t i = 1; i < v.size (); i++)
    v[i].off = v[i-1].off + v[i-1].len;
}

void pack (const atlas::Field & f, size_t lsize, byte * buf)
{
  size_t dsize = f.datatype ().size ();

  switch (f.rank ())
    {
      case 1:
      case 2:
      case 3:
    }
  


}

void gather (const atlas::StructuredGrid & grid, const atlas::grid::Distribution & dist, 
             const atlas::FieldSet & sloc, atlas::FieldSet & sglo)
{
  ATLAS_ASSERT (sloc.size () == sglo.size ());
  size_t nfld = sloc.size ();

  auto & comm = mpi::comm ();
  int nprc = comm.size ();
  int lprc = comm.rank ();

  std::vector<int> count (nprc);
  std::vector<int> ind (grid.size ());

  for (int i = 0; i < grid.size (); i++)
    ind[i] = count[dist.partition (i)]++;

  int max = *std::max_element (count.begin (), count.end ());

  class locprc_t
  {
  public:
    locprc_t () = default;
    locprc_t (const std::initializer_list<int> & l) 
    {
      auto it = l.begin ();
      loc = *it++;
      prc = *it++;
    }
    int loc = std::numeric_limits<int>::min ();
    int prc = std::numeric_limits<int>::min ();
  };

  class locprc2glo_t 
  {
  public:
    locprc2glo_t (int _max, int _nprc) : max (_max), nprc (_nprc)
    {
      data.resize (max * nprc, std::numeric_limits<int>::min ());
    }
    int & operator () (int iprc, int jloc) 
    {
      return data[iprc * max + jloc];
    }
   
    class slice_t
    {
    public:
      slice_t (int _iprc, locprc2glo_t * _locprc2glo) 
             : iprc (_iprc), locprc2glo (_locprc2glo)
      {
      }
      int & operator [] (int jloc)
      {
        return (*locprc2glo) (iprc, jloc);
      }
    private:
      locprc2glo_t * locprc2glo;
      int iprc;
    };

    slice_t operator[] (int iprc)
    {
      return slice_t (iprc, this);
    }
  private:
    int max, nprc;
    std::vector<int> data;
  };

  class glo2locprc_t : public std::vector<locprc_t>
  {
  public:
    glo2locprc_t () = default;
    glo2locprc_t (size_t size) : std::vector<locprc_t> (size)
    {
    }
  };

  class offlen_t
  {
  public:
    size_t off = 0, len = 0;
  };

  class sizmul_t
  {
  public:
    size_t siz = 0, mul = 1;
  };

  locprc2glo_t locprc2glo (max, nprc);
  glo2locprc_t glo2locprc (grid.size ());
  
  for (int i = 0; i < grid.size (); i++)
    {
      locprc2glo[dist.partition (i)][ind[i]] = i;
      glo2locprc[i] = {ind[i], dist.partition (i)};
    }
  
  // Sort fields by owner
  
  std::vector<atlas::Field> floc (nfld);
  std::vector<atlas::Field> fglo (nfld);
  std::vector<int> isort (nfld), owner (nfld);

  for (int jfld = 0; jfld < nfld; jfld++)
    {
      floc[jfld] = sloc[jfld];
      fglo[jfld] = sglo[jfld];
      fglo[jfld].metadata ().get ("owner", owner[jfld]);
    }

  std::iota (std::begin (isort), std::end (isort), 0);
  std::sort (std::begin (isort), std::end (isort), [&owner] (int a, int b) { return owner[a] < owner[b]; });

  owner = reorder (owner, isort);
  floc  = reorder (floc,  isort);
  fglo  = reorder (fglo,  isort);

  // Collect field datatype size & number
  
  std::vector<sizmul_t> fld_info (nfld);

  for (int jfld = 0; jfld < nfld; jfld++)
    {
      auto & f = floc[jfld];
      const auto & ss = f.shape (); 
      for (int i = 1; i < ss.size (); i++)
        fld_info[jfld].mul *= ss[i];
      fld_info[jfld].siz =  f.datatype ().size ();
    }

  // SEND

  std::vector<offlen_t> fld_send (nfld + 1);
  std::vector<offlen_t> prc_send (nprc + 1);

  for (int jfld = 0; jfld < nfld; jfld++)
    {
      // datatype size x leading dimension x product of other dimension
      fld_send[jfld].len = fld_info[jfld].siz * fld_info[jfld].mul;
      prc_send[owner[jfld]].len += fld_send[jfld].len;
    }

  integrate (fld_send);
  integrate (prc_send);

  // Pack send buffer

  std::vector<byte> buf_send (fld_send.back ().off * dist.nb_pts ()[lprc]);

  for (int jfld = 0; jfld < nfld; jfld++)
    pack (floc[jfld], dist.nb_pts ()[lprc], &buf_send[fld_send[jfld].off]);

  // RECV

  std::vector<offlen_t> fld_recv (nfld + 1);

  for (int jfld = 0; jfld < nfld; jfld++)
    if (lprc == owner[jfld])
      fld_recv[jfld].len = fld_info[jfld].siz * fld_info[jfld].mul;

  integrate (fld_recv);

  



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
    atlas::FieldSet sglo, sglo2;

#ifdef UNDEF
    for (int i = 0; i < grid.size (); i++)
      printf (" %8d > %8d, %8d\n", i, prc[i], ind[i]);
#endif


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
        sglo.add (fglo);

        Field fglo2 = Field (name, atlas::array::DataType::kind<T> (), {irank == owner ? grid.size () : 0}); 
        fglo2.metadata ().set ("owner", owner);
        sglo2.add (fglo);

      }

    fs.gather (sloc, sglo);
//  gather (grid, dist, sloc, sglo2);

    for (int i = 0; i < nfields; i++)
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


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
