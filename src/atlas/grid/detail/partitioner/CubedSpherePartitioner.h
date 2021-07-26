/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <vector>

#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class CubedSpherePartitioner : public Partitioner {
public:
    CubedSpherePartitioner();

    CubedSpherePartitioner( int N );  // N is the number of parts (aka MPI tasks)
    CubedSpherePartitioner( int N, const eckit::Parametrisation& );

    CubedSpherePartitioner( const int N, const std::vector<int> & globalProcStartPE,
                                         const std::vector<int> & globalProcEndPE,
                                         const std::vector<int> & nprocx,
                                         const std::vector<int> & nprocy );


    // Node struct that holds the x and y indices (for global, it's longitude and
    // latitude in millidegrees (integers))
    // This structure is used in sorting algorithms, and uses less memory than
    // if x and y were in double precision.
    struct NodeInt {
        int x, y, t;
        int n;
    };

    virtual std::string type() const { return "cubedsphere"; }

private:
    struct CubedSphere {
        idx_t nproc[6];
        idx_t nprocx[6];  // number of PEs in the x direction of xy space on each tile.
        idx_t nprocy[6];  // number of PEs in the y direction of xy space on each tile.
        idx_t globalProcStartPE[6]; // lowest global mpi rank on each tile;
        idx_t globalProcEndPE[6]; // final global mpi rank on each tile;
                   // note that mpi ranks on each tile are vary contiguously from globalProcStartPE to
                   // globalProcEndPE.

        idx_t nx[6], ny[6];  // grid dimensions on each tile - for all cell-centered grids they will be same.

        // the two variables below are for now the main options
        // in the future this will be extended
        idx_t startingCornerOnTile[6] = {0,0,0,0,0,0}  ; // for now bottom left corner (0) default. Could be configurable to
                                       // top left (1), top right(2) bottom right(3)
        idx_t xFirst[6] = {1,1,1,1,1,1}; // if 1 then x is leading index - if 0 y is leading index;
    };

    CubedSphere cubedsphere( const Grid& ) const;

    // Doesn't matter if nodes[] is in degrees or radians, as a sorting
    // algorithm is used internally
    void partition( const CubedSphere& cb, int nb_nodes, NodeInt nodes[], int part[] ) const;

    using Partitioner::partition;
    virtual void partition( const Grid&, int part[] ) const;

    void check() const;

private:
    std::vector<atlas::idx_t> globalProcStartPE_;
    std::vector<atlas::idx_t> globalProcEndPE_;
    std::vector<atlas::idx_t> nprocx_{1,1,1,1,1,1};  // number of ranks in x direction on each tile
    std::vector<atlas::idx_t> nprocy_{1,1,1,1,1,1};  // number of ranks in x direction on each tile
    bool regular_      = true;
    bool cubedsphere_ = true;  // exact (true) or approximate (false) cubedsphere
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
