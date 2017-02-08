#ifndef CheckerBoardPartitioner_h
#define CheckerBoardPartitioner_h

#include <vector>
#include "atlas/grid/partitioners/Partitioner.h"

namespace atlas {
namespace grid { class Grid; }
}

namespace atlas {
namespace grid {
namespace partitioners {

class CheckerBoardPartitioner: public Partitioner
{
public:

  CheckerBoardPartitioner(const grid::Grid&);

  CheckerBoardPartitioner(const grid::Grid&, int N);    // N is the number of parts (aka MPI tasks)

  CheckerBoardPartitioner(const grid::Grid&, int N, int nbands);
  CheckerBoardPartitioner(const grid::Grid&, int N, int nbands, bool checkerboard);

public:

  // Node struct that holds the x and y indices (for global, it's longitude and latitude in millidegrees (integers))
  // This structure is used in sorting algorithms, and uses less memory than
  // if x and y were in double precision.
  struct NodeInt
  {
    int x, y;
    int n;
  };

private:
  // Doesn't matter if nodes[] is in degrees or radians, as a sorting
  // algorithm is used internally
  void partition(int nb_nodes, NodeInt nodes[], int part[]) const;

  void configure_defaults(const grid::Grid&);

private:

  virtual void partition( int part[] ) const;

private:

  int nparts_;  // number of parts
  int nbands_;  // number of bands
  int nx_, ny_;  // grid dimensions
  bool checkerboard_;  // exact (true) or approximate (false) checkerboard

};

} // namespace partitioners
} // namespace grid
} // namespace atlas

#endif // CheckerBoardPartitioner_h
