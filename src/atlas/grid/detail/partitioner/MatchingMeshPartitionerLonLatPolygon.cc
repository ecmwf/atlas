/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerLonLatPolygon.h"

#include <vector>

#include "eckit/config/Resource.h"
#include "eckit/log/ProgressTimer.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/fill.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/PolygonXY.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

namespace {
PartitionerBuilder<MatchingMeshPartitionerLonLatPolygon> __builder("lonlat-polygon");
}

void MatchingMeshPartitionerLonLatPolygon::partition(const Grid& grid, int partitioning[]) const {
    const eckit::mpi::Comm& comm = atlas::mpi::comm();
    const int mpi_rank           = int(comm.rank());
    const int mpi_size           = int(comm.size());

    ATLAS_TRACE("MatchingMeshPartitionerLonLatPolygon::partition");

    ATLAS_ASSERT(grid.domain().global());

    Log::debug() << "MatchingMeshPartitionerLonLatPolygon::partition" << std::endl;

    const util::PolygonXY poly{prePartitionedMesh_.polygon(0)};

    double west = poly.coordinatesMin().x();
    double east = poly.coordinatesMax().x();
    comm.allReduceInPlace(west, eckit::mpi::Operation::MIN);
    comm.allReduceInPlace(east, eckit::mpi::Operation::MAX);

    Projection projection = prePartitionedMesh_.projection();
    omp::fill(partitioning, partitioning + grid.size(), -1);

    auto compute = [&](double west) {
        size_t i = 0;

        for (PointLonLat P : grid.lonlat()) {
            if (partitioning[i] < 0) {
                projection.lonlat2xy(P);
                P.normalise(west);
                partitioning[i] = poly.contains(P) ? mpi_rank : -1;
            }
            ++i;
        }
        // Synchronize partitioning
        comm.allReduceInPlace(partitioning, grid.size(), eckit::mpi::Operation::MAX);

        return *std::min_element(partitioning, partitioning + grid.size());
    };

    int min              = compute(west);
    constexpr double eps = 1.e-10;
    if (min < 0 && east - west > 360. + eps) {
        west = east - 360.;
        min  = compute(west);
    }
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
