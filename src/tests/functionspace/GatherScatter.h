#pragma once


#include "ioFieldDesc.h"
#include "atlas/grid.h"

class GatherScatter
{
private:

  class locprc_t
  {
  public:
    atlas::idx_t loc = std::numeric_limits<atlas::idx_t>::min ();
    atlas::idx_t prc = std::numeric_limits<atlas::idx_t>::min ();
  };

  atlas::idx_t max, nprc;

  std::vector<atlas::gidx_t> _prcloc2glo;
  std::vector<locprc_t> _glo2prcloc;

  const atlas::StructuredGrid & grid;
  const atlas::grid::Distribution & dist;

public:
  GatherScatter (const atlas::StructuredGrid & _grid, const atlas::grid::Distribution & _dist);

  atlas::gidx_t prcloc2glo (atlas::idx_t iprc, atlas::idx_t jloc) const
  {
    return _prcloc2glo[iprc * max + jloc];
  }

  const locprc_t & glo2prcloc (atlas::gidx_t jglo) const
  {
    return _glo2prcloc[jglo];
  }
 
  void gather (std::vector<ioFieldDesc> & floc, std::vector<ioFieldDesc> & fglo) const;
 

};

