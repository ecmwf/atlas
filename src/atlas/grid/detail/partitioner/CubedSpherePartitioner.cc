/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "CubedSpherePartitioner.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <functional>
#include <iostream>
#include <vector>

#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

namespace {

bool isNearInt(double value) {
    const double diff = value - floor(value);
    return (diff <= std::numeric_limits<double>::epsilon() || diff >= (1.0 - std::numeric_limits<double>::epsilon()));
}

bool isConfigSufficient(const eckit::Parametrisation& config) {
    return !(config.has("starting rank on tile")) || !(config.has("final rank on tile")) || !(config.has("nprocx")) ||
           !(config.has("nprocy"));
}

std::vector<std::vector<atlas::idx_t>> createOffset(const std::array<std::size_t, 6>& ngpt,
                                                    const std::array<idx_t, 6> nproc1D,
                                                    const std::array<idx_t, 6> maxDim1D,
                                                    const std::array<idx_t, 6> maxDim2D) {
    std::vector<std::vector<atlas::idx_t>> offset;

    idx_t totalNproc1D = std::accumulate(std::begin(nproc1D), std::end(nproc1D), 0);

    // ngpbt - number of grid points per band and tile.
    std::vector<std::size_t> ngpbt(static_cast<std::size_t>(totalNproc1D), 0);

    idx_t j(0);
    for (std::size_t t = 0; t < 6; ++t) {
        idx_t jt          = static_cast<idx_t>(j);
        std::size_t nproc = static_cast<std::size_t>(nproc1D[t]);
        std::vector<idx_t> offsetPerTile(nproc + 1, 0);

        for (std::size_t proc1D = 0; proc1D < nproc; ++proc1D, ++j) {
            ngpbt[static_cast<size_t>(j)] =
                proc1D < nproc - 1
                    ? ngpt[t] / nproc
                    : ngpt[t] - static_cast<std::size_t>(std::accumulate(ngpbt.begin() + jt, ngpbt.begin() + j, 0));

            offsetPerTile[proc1D] =
                proc1D == 0 ? 0 : std::accumulate(ngpbt.begin() + jt, ngpbt.begin() + j, 0) / maxDim2D[t];
        }
        offsetPerTile[nproc] = maxDim1D[t];
        offset.push_back(offsetPerTile);
    }

    return offset;
}

}  // namespace

CubedSpherePartitioner::CubedSpherePartitioner(): Partitioner() {}

CubedSpherePartitioner::CubedSpherePartitioner(int N): Partitioner(N, util::NoConfig()), regular_{true} {}

CubedSpherePartitioner::CubedSpherePartitioner(const eckit::Parametrisation& config):
    Partitioner(config), regular_{isConfigSufficient(config)} {
    if (config.has("starting rank on tile") && config.has("final rank on tile") && config.has("nprocx") &&
        config.has("nprocy")) {
        config.get("starting rank on tile", globalProcStartPE_);
        config.get("final rank on tile", globalProcEndPE_);
        config.get("nprocx", nprocx_);
        config.get("nprocy", nprocy_);
    }
}

CubedSpherePartitioner::CubedSpherePartitioner(int N, const eckit::Parametrisation& config):
    Partitioner(N, config), regular_{isConfigSufficient(config)} {
    if (config.has("starting rank on tile") && config.has("final rank on tile") && config.has("nprocx") &&
        config.has("nprocy")) {
        config.get("starting rank on tile", globalProcStartPE_);
        config.get("final rank on tile", globalProcEndPE_);
        config.get("nprocx", nprocx_);
        config.get("nprocy", nprocy_);
    }
}

CubedSpherePartitioner::CubedSpherePartitioner(const int N, const std::vector<int>& globalProcStartPE,
                                               const std::vector<int>& globalProcEndPE, const std::vector<int>& nprocx,
                                               const std::vector<int>& nprocy):
    Partitioner(N, util::NoConfig()),
    globalProcStartPE_(globalProcStartPE.begin(), globalProcStartPE.end()),
    globalProcEndPE_(globalProcEndPE.begin(), globalProcEndPE.end()),
    nprocx_{nprocx.begin(), nprocx.end()},
    nprocy_{nprocy.begin(), nprocy.end()},
    regular_{false} {}


CubedSpherePartitioner::CubedSphere CubedSpherePartitioner::cubedsphere(const Grid& grid) const {
    // grid dimensions
    const CubedSphereGrid cg(grid);
    if (!cg) {
        throw_Exception("CubedSphere Partitioner only works for cubed sphere grids.", Here());
    }

    CubedSphere cb;

    for (std::size_t t = 0; t < 6; ++t) {
        cb.nx[t] = cg.N();
        cb.ny[t] = cg.N();
    }

    atlas::idx_t nparts = nb_partitions();

    if (regular_) {
        // share PEs around tiles
        // minRanksPerTile

        idx_t ranksPerTile = nparts / 6;

        idx_t reminder = nparts - 6 * ranksPerTile;

        for (std::size_t t = 0; t < 6; ++t) {
            cb.nproc[t] = ranksPerTile;
        }

        // round-robin;
        std::size_t t{0};
        while (reminder > 0) {
            if (t == 6) {
                t = 0;
            }
            cb.nproc[t] += 1;
            t += 1;
            reminder -= 1;
        }

        // now need to specify nprocx and nprocy.
        // nproc is 0 for the tile we use default nprocx and nprocy = 1

        // if we can square-root nproc and get an integer
        // we use that for nprocx and nprocy
        // otherwise we split just in nprocx and keep nprocy =1;

        for (std::size_t t = 0; t < 6; ++t) {
            if (cb.nproc[t] > 0) {
                double sq = std::sqrt(static_cast<double>(cb.nproc[t]));
                if (isNearInt(sq)) {
                    cb.nprocx[t] = static_cast<idx_t>(std::round(sq));
                    cb.nprocy[t] = static_cast<idx_t>(std::round(sq));
                }
                else {
                    cb.nprocx[t] = 1;
                    cb.nprocy[t] = cb.nproc[t];
                }
            }
        }
        cb.globalProcStartPE[0] = 0;
        cb.globalProcEndPE[0]   = cb.nproc[0] - 1;

        for (size_t t = 1; t < 6; ++t) {
            if (cb.nproc[t] == 0) {
                cb.globalProcStartPE[t] = cb.globalProcEndPE[t - 1];
                cb.globalProcEndPE[t]   = cb.globalProcEndPE[t - 1];
            }
            else {
                cb.globalProcStartPE[t] = cb.globalProcEndPE[t - 1] + 1;
                cb.globalProcEndPE[t]   = cb.globalProcStartPE[t] + cb.nproc[t] - 1;
            }
        }
    }
    else {
        for (size_t t = 0; t < 6; ++t) {
            cb.globalProcStartPE[t] = globalProcStartPE_[t];
            cb.globalProcEndPE[t]   = globalProcEndPE_[t];
            cb.nprocx[t]            = nprocx_[t];
            cb.nprocy[t]            = nprocy_[t];
        }
        cb.nproc[0] = globalProcEndPE_[0] + 1;
        for (std::size_t t = 1; t < 6; ++t) {
            cb.nproc[t] = cb.globalProcStartPE[t] == cb.globalProcEndPE[t - 1]
                              ? cb.nproc[t] = 0
                              : cb.globalProcEndPE[t] - cb.globalProcStartPE[t] + 1;
        }
    }
    return cb;
}

void CubedSpherePartitioner::partition(CubedSphere& cb, const int nb_cells, const CellInt cells[], int part[]) const {
    // ngpt number of grid points per tile (this is for cell centers where the number of cells are the same
    std::size_t tileCells = static_cast<std::size_t>(nb_cells / 6);
    std::array<std::size_t, 6> ngpt{tileCells, tileCells, tileCells, tileCells, tileCells, tileCells};

    // xoffset - the xoffset per band per tile (including cb.nx[t])
    //std::vector<std::vector<atlas::idx_t>> xoffset(createOffset(ngpt, cb.nprocx, cb.nx, cb.ny));

    cb.xoffset = createOffset(ngpt, cb.nprocx, cb.nx, cb.ny);

    // yoffset - the yoffset per band per tile (including cb.ny[t])
    cb.yoffset = createOffset(ngpt, cb.nprocy, cb.ny, cb.nx);

    int p;  // partition index
    // loop over all data tile by tile
    for (int n = 0; n < nb_cells; ++n) {
        std::size_t t = static_cast<std::size_t>(cells[n].t);
        p             = cb.globalProcStartPE[t];
        for (std::size_t yproc = 0; yproc < static_cast<std::size_t>(cb.nprocy[t]); ++yproc) {
            for (std::size_t xproc = 0; xproc < static_cast<std::size_t>(cb.nprocx[t]); ++xproc) {
                if ((cells[n].y >= cb.yoffset[t][yproc]) && (cells[n].y < cb.yoffset[t][yproc + 1]) &&
                    (cells[n].x >= cb.xoffset[t][xproc]) && (cells[n].x < cb.xoffset[t][xproc + 1])) {
                    part[n] = p;
                }
                ++p;
            }
        }
    }
    ATLAS_ASSERT(part[nb_cells - 1] == nb_partitions() - 1, "number of partitions created not what is expected");
}

void CubedSpherePartitioner::partition(const Grid& grid, int part[]) const {
    if (nb_partitions() == 1)  // trivial solution, so much faster
    {
        for (idx_t j = 0; j < grid.size(); ++j) {
            part[j] = 0;
        }
    }
    else {
        auto cb = cubedsphere(grid);

        std::vector<CellInt> cells(static_cast<std::size_t>(grid.size()));
        std::size_t n(0);

        for (std::size_t it = 0; it < 6; ++it) {
            for (idx_t iy = 0; iy < cb.ny[it]; ++iy) {
                for (idx_t ix = 0; ix < cb.nx[it]; ++ix) {
                    cells[n].t = static_cast<int>(it);
                    cells[n].x = static_cast<int>(ix);
                    cells[n].y = static_cast<int>(iy);
                    cells[n].n = static_cast<int>(n);
                    ++n;
                }
            }
        }
        partition(cb, grid.size(), cells.data(), part);
    }
}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas

namespace {
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::CubedSpherePartitioner>
    __CubedSphere("cubedsphere");
}
