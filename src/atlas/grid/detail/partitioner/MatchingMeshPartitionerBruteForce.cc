/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerBruteForce.h"

#include <vector>

#include "eckit/log/ProgressTimer.h"

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/grid.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

namespace {

PartitionerBuilder<MatchingMeshPartitionerBruteForce> __builder("brute-force");

double dot_sign(const double& Ax, const double& Ay, const double& Bx, const double& By, const double& Cx,
                const double& Cy) {
    return (Ax - Cx) * (By - Cy) - (Bx - Cx) * (Ay - Cy);
}

/// Point-in-triangle test
/// @note "dot product" test, equivalent to half-plane check
/// @see
/// http://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
bool point_in_triangle(const PointLonLat& P, const PointLonLat& A, const PointLonLat& B, const PointLonLat& C) {
    // Compute signs of dot products
    const bool b1 = dot_sign(P.lon(), P.lat(), A.lon(), A.lat(), B.lon(), B.lat()) < 0,
               b2 = dot_sign(P.lon(), P.lat(), B.lon(), B.lat(), C.lon(), C.lat()) < 0,
               b3 = dot_sign(P.lon(), P.lat(), C.lon(), C.lat(), A.lon(), A.lat()) < 0;

    // Check if point is in triangle
    // - check for b1 && b2 && b3 only works for triangles ordered
    // counter-clockwise
    // - check for b1==b2 && b2==b3 works for both, equivalent to "all points on
    // the same side"
    return (b1 == b2) && (b2 == b3);
}

/// Point-in-quadrilateral test
/// @note limited to convex quadrilaterals, @see
/// https://en.wikipedia.org/wiki/Convex_polygon
bool point_in_quadrilateral(const PointLonLat& P, const PointLonLat& A, const PointLonLat& B, const PointLonLat& C,
                            const PointLonLat& D) {
    const bool b1 = dot_sign(P.lon(), P.lat(), A.lon(), A.lat(), B.lon(), B.lat()) < 0,
               b2 = dot_sign(P.lon(), P.lat(), B.lon(), B.lat(), C.lon(), C.lat()) < 0,
               b3 = dot_sign(P.lon(), P.lat(), C.lon(), C.lat(), D.lon(), D.lat()) < 0,
               b4 = dot_sign(P.lon(), P.lat(), D.lon(), D.lat(), A.lon(), A.lat()) < 0;
    return (b1 == b2) && (b2 == b3) && (b3 == b4);
}

}  // namespace

void MatchingMeshPartitionerBruteForce::partition(const Grid& grid, int partitioning[]) const {
    ATLAS_TRACE("MatchingMeshPartitionerBruteForce::partition");

    const auto& comm   = mpi::comm(prePartitionedMesh_.mpi_comm());
    const int mpi_rank = int(comm.rank());
    const int mpi_size = int(comm.size());

    // Point coordinates
    // - use a bounding box to quickly discard points,
    // - except when that is above/below bounding box but poles should be included
    //   FIXME: THIS IS A HACK! the coordinates include North/South Pole
    //   (first/last partitions only)

    ATLAS_ASSERT(grid.domain().global());
    bool includesNorthPole = (mpi_rank == 0);
    bool includesSouthPole = (mpi_rank == (int(comm.size()) - 1));

    ATLAS_ASSERT(prePartitionedMesh_.nodes().size());
    auto lonlat_src = array::make_view<double, 2>(prePartitionedMesh_.nodes().lonlat());

    std::vector<PointLonLat> coordinates;
    coordinates.reserve(lonlat_src.shape(0));
    PointLonLat coordinatesMin = PointLonLat(lonlat_src(0, LON), lonlat_src(0, LAT));
    PointLonLat coordinatesMax = coordinatesMin;

    for (idx_t i = 0; i < lonlat_src.shape(0); ++i) {
        PointLonLat A(lonlat_src(i, LON), lonlat_src(i, LAT));
        coordinatesMin = PointLonLat::componentsMin(coordinatesMin, A);
        coordinatesMax = PointLonLat::componentsMax(coordinatesMax, A);
        coordinates.emplace_back(A);
    }

    {
        eckit::ProgressTimer timer("Partitioning target", grid.size(), "point", double(10), atlas::Log::trace());
        auto grid_iter = grid.lonlat().begin();
        for (idx_t i = 0; i < grid.size(); ++i, ++grid_iter) {
            partitioning[i] = -1;
            const PointLonLat& P = *grid_iter;

            if (coordinatesMin[LON] <= P[LON] && P[LON] <= coordinatesMax[LON] && coordinatesMin[LAT] <= P[LAT] &&
                P[LAT] <= coordinatesMax[LAT]) {
                const mesh::Cells& elements_src = prePartitionedMesh_.cells();
                const idx_t nb_types            = elements_src.nb_types();
                bool found                      = false;

                for (idx_t t = 0; t < nb_types && !found; ++t) {
                    idx_t idx[4];
                    const mesh::Elements& elements      = elements_src.elements(t);
                    const mesh::BlockConnectivity& conn = elements.node_connectivity();

                    const idx_t nb_nodes = elements.nb_nodes();
                    ATLAS_ASSERT((nb_nodes == 3 && elements.name() == "Triangle") ||
                                 (nb_nodes == 4 && elements.name() == "Quadrilateral"));

                    for (idx_t j = 0; j < elements.size() && !found; ++j) {
                        idx[0] = conn(j, 0);
                        idx[1] = conn(j, 1);
                        idx[2] = conn(j, 2);
                        idx[3] = nb_nodes > 3 ? conn(j, 3) : 0;

                        if ((elements.name() == "Triangle" &&
                             point_in_triangle(P, coordinates[idx[0]], coordinates[idx[1]], coordinates[idx[2]])) ||
                            (elements.name() == "Quadrilateral" &&
                             point_in_quadrilateral(P, coordinates[idx[0]], coordinates[idx[1]], coordinates[idx[2]],
                                                    coordinates[idx[3]]))) {
                            partitioning[i] = mpi_rank;
                            found           = true;
                        }
                    }
                }
            }
            else if ((includesNorthPole && P[LAT] > coordinatesMax[LAT]) ||
                     (includesSouthPole && P[LAT] < coordinatesMin[LAT])) {
                partitioning[i] = mpi_rank;
            }
        }
    }

    // Synchronize the partitioning and return a grid partitioner
    comm.allReduceInPlace(partitioning, grid.size(), eckit::mpi::Operation::MAX);
    const int min = *std::min_element(partitioning, partitioning + grid.size());
    if (min < 0) {
        throw_Exception(
            "Could not find partition for target node (source "
            "mesh does not contain all target grid points)",
            Here());
    }
}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
