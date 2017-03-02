#pragma once

#include <vector>
#include "atlas/grid/detail/partitioners/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioners {

class CheckerBoardPartitioner: public Partitioner {

public:

    CheckerBoardPartitioner(const Grid&);

    CheckerBoardPartitioner(const Grid&, int N);    // N is the number of parts (aka MPI tasks)

    CheckerBoardPartitioner(const Grid&, int N, int nbands);
    CheckerBoardPartitioner(const Grid&, int N, int nbands, bool checkerboard);

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

    void configure_defaults(const Grid&);

    virtual void partition( int part[] ) const;

private:

    size_t nparts_;  // number of parts
    size_t nbands_;  // number of bands
    size_t nx_, ny_;  // grid dimensions
    bool checkerboard_;  // exact (true) or approximate (false) checkerboard

};

} // namespace partitioners
} // namespace detail
} // namespace grid
} // namespace atlas
