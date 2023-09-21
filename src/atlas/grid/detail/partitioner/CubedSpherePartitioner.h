/*
 * (C) Crown Copyright Met Office 2021
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include <vector>

#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class CubedSpherePartitioner : public Partitioner {
public:
    CubedSpherePartitioner();

    CubedSpherePartitioner(int N);  // N is the number of parts (aka MPI tasks)
    CubedSpherePartitioner(int N, const eckit::Parametrisation&);
    CubedSpherePartitioner(const eckit::Parametrisation&);

    CubedSpherePartitioner(const int N, const std::vector<int>& globalProcStartPE,
                           const std::vector<int>& globalProcEndPE, const std::vector<int>& nprocx,
                           const std::vector<int>& nprocy);

    CubedSpherePartitioner(const int N, const bool regularGrid);


    // Cell struct that holds the x and y and t indices
    // This could be replaced using the tij iterator that is part of the grid class
    struct CellInt {
        int x, y, t;
        int n;
    };

    struct CubedSphere {
        std::array<atlas::idx_t, 6> nproc;
        std::array<atlas::idx_t, 6> nprocx;             // number of PEs in the x direction of xy space on each tile.
        std::array<atlas::idx_t, 6> nprocy;             // number of PEs in the y direction of xy space on each tile.
        std::array<atlas::idx_t, 6> globalProcStartPE;  // lowest global mpi rank on each tile;
        std::array<atlas::idx_t, 6> globalProcEndPE;    // final global mpi rank on each tile;
            // note that mpi ranks on each tile are vary contiguously from globalProcStartPE to
            // globalProcEndPE.

        // grid dimensions on each tile - for all cell-centered grids they will be same.
        std::array<atlas::idx_t, 6> nx;
        std::array<atlas::idx_t, 6> ny;

        // these are the offsets in the x and y directions
        // they are allocated in "void partition(CubedSphere& cb, int nb_nodes, CellInt nodes[], int part[] );"
        std::vector<std::vector<atlas::idx_t>> xoffset;
        std::vector<std::vector<atlas::idx_t>> yoffset;

        // the two variables below are for now the main options
        // in the future this will be extended
        std::array<atlas::idx_t, 6> startingCornerOnTile;
        // for now bottom left corner (0) default.
        // Could be configurable to top left (1), top right(2) bottom right(3)

        std::array<atlas::idx_t, 6> xFirst;  // if 1 then x is leading index - if 0 y is leading index;

        CubedSphere() {
            nprocx               = std::array<atlas::idx_t, 6>{1, 1, 1, 1, 1, 1};
            nprocy               = std::array<atlas::idx_t, 6>{1, 1, 1, 1, 1, 1};
            startingCornerOnTile = std::array<atlas::idx_t, 6>{0, 0, 0, 0, 0, 0};
            xFirst               = std::array<atlas::idx_t, 6>{1, 1, 1, 1, 1, 1};
        }
    };

    CubedSphere cubedsphere(const Grid&) const;

    void partition(CubedSphere& cb, const int nb_nodes, const CellInt nodes[], int part[]) const;

    virtual std::string type() const { return "cubedsphere"; }

private:
    using Partitioner::partition;
    virtual void partition(const Grid&, int part[]) const;

    void check() const;

private:
    std::vector<atlas::idx_t> globalProcStartPE_{0, 0, 0, 0, 0, 0};
    std::vector<atlas::idx_t> globalProcEndPE_{0, 0, 0, 0, 0, 0};
    std::vector<atlas::idx_t> nprocx_{1, 1, 1, 1, 1, 1};  // number of ranks in x direction on each tile
    std::vector<atlas::idx_t> nprocy_{1, 1, 1, 1, 1, 1};  // number of ranks in x direction on each tile
    bool regular_ = true;                                 // regular algorithm for partitioning.
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
