#pragma once


#include "ioFieldDesc.h"
#include "atlas/grid.h"

class GatherScatter
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
  GatherScatter (const atlas::StructuredGrid & _grid, const atlas::grid::Distribution & _dist);

  int prcloc2glo (int iprc, int jloc) const
  {
    return _prcloc2glo[iprc * max + jloc];
  }

  const locprc_t & glo2prcloc (int jglo) const
  {
    return _glo2prcloc[jglo];
  }
 
  void gather (std::vector<ioFieldDesc> & floc, std::vector<ioFieldDesc> & fglo) const;
 

};

