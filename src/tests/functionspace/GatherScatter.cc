#include "GatherScatter.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"

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

template <typename I, typename F>
std::vector<atlas::idx_t> grep (I n, F f)
{
  std::vector<atlas::idx_t> r (n), g;
  std::iota (r.begin (), r.end (), 0);
  std::copy_if (r.begin (), r.end (), std::back_inserter (g), f);
  return g;
}

};

GatherScatter::GatherScatter (const atlas::StructuredGrid & _grid, const atlas::grid::Distribution & _dist)
 : grid (_grid), dist (_dist)
{
ATLAS_TRACE_SCOPE ("GatherScatter::GatherScatter")
{

  nprc = dist.nb_partitions ();

  std::vector<atlas::gidx_t> count (nprc);
  std::vector<atlas::idx_t> ind (grid.size ());

  for (atlas::gidx_t i = 0; i < grid.size (); i++)
    ind[i] = count[dist.partition (i)]++;

  max = *std::max_element (count.begin (), count.end ());

  _prcloc2glo.resize (max * nprc, std::numeric_limits<atlas::gidx_t>::min ());
  _glo2prcloc.resize (grid.size ());

  for (atlas::gidx_t i = 0; i < grid.size (); i++)
    {
      _prcloc2glo[max * dist.partition (i) + ind[i]] = i;
      _glo2prcloc[i].loc = ind[i];
      _glo2prcloc[i].prc = dist.partition (i);
    }

}
}

void GatherScatter::gather (ioFieldDesc_v & floc, ioFieldDesc_v & fglo) const
{
ATLAS_TRACE_SCOPE ("GatherScatter::gather")
{

  ATLAS_ASSERT (floc.size () == fglo.size ());
  atlas::idx_t nfld = floc.size ();

  auto & comm = eckit::mpi::comm ();
  atlas::idx_t nprc = comm.size ();
  atlas::idx_t lprc = comm.rank ();

  atlas::idx_t ldim = dist.nb_pts ()[lprc];

  // Sort fields by owner
  
  std::vector<atlas::idx_t> isort (nfld);

  std::iota (std::begin (isort), std::end (isort), 0);
  std::sort (std::begin (isort), std::end (isort), 
             [&fglo] (atlas::idx_t a, atlas::idx_t b) 
             { return fglo[a].owner () < fglo[b].owner (); });

  floc  = reorder (floc,  isort);
  fglo  = reorder (fglo,  isort);

  // SEND

  offlen_v fld_loc (nfld + 1);
  offlen_v prc_loc (nprc + 1);

  for (atlas::idx_t jfld = 0; jfld < nfld; jfld++)
    {
      atlas::idx_t owner = fglo[jfld].owner ();
      fld_loc[jfld].len = floc[jfld].size ();
      prc_loc[owner].len += floc[jfld].size ();
    }

  integrate (fld_loc);
  integrate (prc_loc);

  // RECV

  offlen_v prc_glo (nprc + 1);
  offlen_v fld_glo (nfld + 1);

  for (atlas::idx_t jfld = 0; jfld < nfld; jfld++)
    if (lprc == fglo[jfld].owner ())
      fld_glo[jfld].len = fglo[jfld].dlen ();

  integrate (fld_glo);

  for (atlas::idx_t iprc = 0; iprc < nprc; iprc++)
     prc_glo[iprc].len = dist.nb_pts ()[iprc] * fld_glo.back ().off;

  integrate (prc_glo);

  // Pack send buffer

  std::vector<byte> buf_loc (fld_loc.back ().off);

  ATLAS_TRACE_SCOPE ("Pack")
  {
#pragma omp parallel for
    for (atlas::idx_t jfld = 0; jfld < nfld; jfld++)
      {
        auto & f = floc[jfld];
        byte * buffer = &buf_loc[fld_loc[jfld].off];
        for (atlas::idx_t i = 0; i < f.ldim (); i++)
        for (int j = 0; j < f.dlen (); j++)
          buffer[i*f.dlen ()+j] = f (i, j);
      }
  }

  std::vector<byte> buf_glo (prc_glo.back ().off);

  ATLAS_TRACE_SCOPE ("SEND/RECV")
  {
    std::vector<eckit::mpi::Request> rqr;
   
    for (atlas::idx_t iprc = 0; iprc < nprc; iprc++)
      if (prc_glo[iprc].len > 0)
        rqr.push_back (comm.iReceive (&buf_glo[prc_glo[iprc].off], 
                                      prc_glo[iprc].len, iprc, 100));
   
    comm.barrier ();
   
    std::vector<eckit::mpi::Request> rqs;
   
    for (atlas::idx_t iprc = 0; iprc < nprc; iprc++)
      if (prc_loc[iprc].len > 0)
        rqs.push_back (comm.iSend (&buf_loc[prc_loc[iprc].off], 
                                   prc_loc[iprc].len, iprc, 100));
   
    for (auto & r : rqr)
      comm.wait (r);
   
    for (auto & r : rqs)
      comm.wait (r);
  }

  ATLAS_TRACE_SCOPE ("Unpack")
  {
    std::vector<atlas::idx_t> prcs = grep (nprc, 
       [&prc_glo] (atlas::idx_t i) { return prc_glo[i].len > 0; });
    std::vector<atlas::idx_t> flds = grep (nfld, 
       [&fld_glo] (atlas::idx_t i) { return fld_glo[i].len > 0; });

#ifdef UNDEF
    for (atlas::idx_t iprc = 0; iprc < nprc; iprc++)
      if (prc_glo[iprc].len > 0)
        {
          for (atlas::idx_t jfld = 0; jfld < nfld; jfld++)
            {
              if (fld_glo[jfld].len > 0)
#endif

    for (int ii = 0; ii < prcs.size (); ii++)
    for (int jj = 0; jj < flds.size (); jj++)
                {
                  const atlas::idx_t iprc = prcs[ii];
                  const atlas::idx_t jfld = flds[jj];
                  const atlas::idx_t ngptot = dist.nb_pts ()[iprc];
                  const size_t off = prc_glo[iprc].off + ngptot * fld_glo[jfld].off;
                  auto & f = fglo[jfld];
                  const size_t len = fld_glo[jfld].len;
#pragma omp parallel for 
                  for (atlas::idx_t jloc = 0; jloc < ngptot; jloc++)
                  for (int j = 0; j < len; j++)
                    {
                      size_t k = jloc * fld_glo[jfld].len + j;
                      f (prcloc2glo (iprc, jloc), j) = buf_glo[off+k];
                    }

                }

  }


}
}

