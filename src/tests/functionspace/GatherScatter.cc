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

void GatherScatter::reOrderFields (ioFieldDesc_v & floc, ioFieldDesc_v & fglo) const
{
  ATLAS_ASSERT (floc.size () == fglo.size ());
  atlas::idx_t nfld = fglo.size ();  

  // Sort fields by owner
  
  std::vector<atlas::idx_t> isort (nfld);

  std::iota (std::begin (isort), std::end (isort), 0);
  std::sort (std::begin (isort), std::end (isort), 
             [&fglo] (atlas::idx_t a, atlas::idx_t b) 
             { return fglo[a].owner () < fglo[b].owner (); });

  floc  = reorder (floc,  isort);
  fglo  = reorder (fglo,  isort);

  for (int jfld = 0; jfld < nfld; jfld++)
    floc[jfld].owner () = fglo[jfld].owner ();

}

void GatherScatter::computeTLoc (const ioFieldDesc_v & floc, fldprc_t & tloc) const
{
  auto & comm = eckit::mpi::comm ();
  atlas::idx_t nprc = comm.size ();
  atlas::idx_t nfld = floc.size ();

  tloc.fld.resize (nfld + 1);
  tloc.prc.resize (nprc + 1);

  for (atlas::idx_t jfld = 0; jfld < nfld; jfld++)
    {
      atlas::idx_t owner = floc[jfld].owner ();
      tloc.fld[jfld].len = floc[jfld].size ();
      tloc.prc[owner].len += floc[jfld].size ();
    }

  integrate (tloc.fld);
  integrate (tloc.prc);
}

void GatherScatter::computeTGlo (const ioFieldDesc_v & fglo, fldprc_t & tglo) const
{
  auto & comm = eckit::mpi::comm ();
  atlas::idx_t nprc = comm.size ();
  atlas::idx_t lprc = comm.rank ();
  atlas::idx_t nfld = fglo.size ();

  tglo.prc.resize (nprc + 1);
  tglo.fld.resize (nfld + 1);

  for (atlas::idx_t jfld = 0; jfld < nfld; jfld++)
    if (lprc == fglo[jfld].owner ())
      tglo.fld[jfld].len = fglo[jfld].dlen ();

  integrate (tglo.fld);

  for (atlas::idx_t iprc = 0; iprc < nprc; iprc++)
     tglo.prc[iprc].len = dist.nb_pts ()[iprc] * tglo.fld.back ().off;

  integrate (tglo.prc);
}


template <typename A>
void GatherScatter::processLocBuffer (ioFieldDesc_v & floc, const fldprc_t & tloc,
                                      byte_v & buf_loc, A a) const
{
  atlas::idx_t nfld = floc.size ();

  ATLAS_TRACE_SCOPE ("GatherScatter::processLocBuffer")
  {
#pragma omp parallel for 
    for (atlas::idx_t jfld = 0; jfld < nfld; jfld++)
      {
        auto & f = floc[jfld];
        byte * buffer = &buf_loc[tloc.fld[jfld].off];
        const size_t dlen = f.dlen ();
        const size_t ldim = f.ldim ();
        for (atlas::idx_t i = 0; i < ldim; i++)
        for (int j = 0; j < dlen; j++)
          a (buffer[i*dlen+j], f (i, j));
      }
  }

}

template <typename A>
void GatherScatter::processGloBuffer (ioFieldDesc_v & fglo, const fldprc_t & tglo,
                                      byte_v & buf_glo, A a) const
{
  auto & comm = eckit::mpi::comm ();
  atlas::idx_t nfld = fglo.size ();
  atlas::idx_t nprc = comm.size ();
  
  ATLAS_TRACE_SCOPE ("GatherScatter::processGloBuffer")
  {
    std::vector<atlas::idx_t> prcs = grep (nprc, 
       [&tglo] (atlas::idx_t i) { return tglo.prc[i].len > 0; });
    std::vector<atlas::idx_t> flds = grep (nfld, 
       [&tglo] (atlas::idx_t i) { return tglo.fld[i].len > 0; });

#ifdef UNDEF
    for (atlas::idx_t iprc = 0; iprc < nprc; iprc++)
      if (tglo.prc[iprc].len > 0)
        {
          for (atlas::idx_t jfld = 0; jfld < nfld; jfld++)
            {
              if (tglo.fld[jfld].len > 0)
#endif

    for (int ii = 0; ii < prcs.size (); ii++)
    for (int jj = 0; jj < flds.size (); jj++)
                {
                  const atlas::idx_t iprc = prcs[ii];
                  const atlas::idx_t jfld = flds[jj];
                  const atlas::idx_t ngptot = dist.nb_pts ()[iprc];
                  const size_t off = tglo.prc[iprc].off + ngptot * tglo.fld[jfld].off;
                  auto & f = fglo[jfld];
                  const size_t len = tglo.fld[jfld].len;
#pragma omp parallel for 
                  for (atlas::idx_t jloc = 0; jloc < ngptot; jloc++)
                  for (int j = 0; j < len; j++)
                    {
                      size_t k = jloc * tglo.fld[jfld].len + j;
                      a (buf_glo[off+k], f (prcloc2glo (iprc, jloc), j));
                    }

                }

  }
}

void GatherScatter::gather (const ioFieldDesc_v & _floc, ioFieldDesc_v & fglo) const
{

ioFieldDesc_v floc = _floc;

ATLAS_TRACE_SCOPE ("GatherScatter::gather")
{
  auto & comm = eckit::mpi::comm ();
  atlas::idx_t nfld = floc.size ();
  atlas::idx_t nprc = comm.size ();
  atlas::idx_t lprc = comm.rank ();
  atlas::idx_t ldim = dist.nb_pts ()[lprc];

  reOrderFields (floc, fglo);

  fldprc_t tloc, tglo;

  computeTLoc (floc, tloc);
  computeTGlo (fglo, tglo);

  byte_v buf_loc (tloc.fld.back ().off);

  processLocBuffer (floc, tloc, buf_loc, [] (byte & a, const byte & b) { a = b; });

  byte_v buf_glo (tglo.prc.back ().off);

  ATLAS_TRACE_SCOPE ("SEND/RECV")
  {
    std::vector<eckit::mpi::Request> rqr;
   
    for (atlas::idx_t iprc = 0; iprc < nprc; iprc++)
      if (tglo.prc[iprc].len > 0)
        rqr.push_back (comm.iReceive (&buf_glo[tglo.prc[iprc].off], 
                                      tglo.prc[iprc].len, iprc, 100));
   
    comm.barrier ();
   
    std::vector<eckit::mpi::Request> rqs;
   
    for (atlas::idx_t iprc = 0; iprc < nprc; iprc++)
      if (tloc.prc[iprc].len > 0)
        rqs.push_back (comm.iSend (&buf_loc[tloc.prc[iprc].off], 
                                   tloc.prc[iprc].len, iprc, 100));
   
    for (auto & r : rqr)
      comm.wait (r);
   
    for (auto & r : rqs)
      comm.wait (r);
  }

  processGloBuffer (fglo, tglo, buf_glo, [] (const byte & a, byte & b) { b = a; });

}
}

