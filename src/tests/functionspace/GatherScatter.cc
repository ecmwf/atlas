#include "GatherScatter.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"

namespace
{

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

};

GatherScatter::GatherScatter (const atlas::StructuredGrid & _grid, const atlas::grid::Distribution & _dist)
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

void GatherScatter::gather (std::vector<ioFieldDesc> & floc, std::vector<ioFieldDesc> & fglo) const
{
  ATLAS_ASSERT (floc.size () == fglo.size ());
  size_t nfld = floc.size ();

  auto & comm = eckit::mpi::comm ();
  int nprc = comm.size ();
  int lprc = comm.rank ();

  class offlen_t
  {
  public:
    size_t off = 0, len = 0;
  };

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

  integrate (fld_send);
  integrate (prc_send);

  // Pack send buffer

  std::vector<byte> buf_send (fld_send.back ().off);

  for (int jfld = 0; jfld < nfld; jfld++)
    {
      auto & f = floc[jfld];
      byte * buffer = &buf_send[fld_send[jfld].off];
      for (int i = 0; i < f.ldim (); i++)
      for (int j = 0; j < f.dlen (); j++)
        buffer[i*f.dlen ()+j] = f (i, j);
    }

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
            auto & f = fglo[jfld];
            if (fld_recv[jfld].len > 0)
              {
                for (int jloc = 0, k = 0; jloc < dist.nb_pts ()[iprc]; jloc++)
                for (int j = 0; j < fld_recv[jfld].len; j++, k++)
                  f (prcloc2glo (iprc, jloc), j) = buf_recv[off+k];
                off += dist.nb_pts ()[iprc] * f.dlen ();
              }
          }
      }
}

