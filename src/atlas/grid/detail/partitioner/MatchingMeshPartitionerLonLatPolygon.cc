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

#include <iomanip>
#include <sstream>
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

    int min              = compute(east - 360.);
    constexpr double eps = 1.e-10;
    bool second_try      = [&]() {
        if (min < 0 && east - west > 360. + eps) {
            min = compute(west - eps);
            return true;
        }
        return false;
    }();
    if (min < 0) {
        size_t i            = 0;
        size_t max_failures = grid.size();
        std::vector<size_t> failed_index;
        std::vector<PointLonLat> failed_lonlat;
        failed_index.reserve(max_failures);
        failed_lonlat.reserve(max_failures);
        for (PointLonLat P : grid.lonlat()) {
            if (partitioning[i] < 0) {
                failed_index.emplace_back(i);
                failed_lonlat.emplace_back(P);
            }
            ++i;
        }
        size_t nb_failures = failed_index.size();
        std::stringstream err;
        err.precision(20);
        err << "Could not find partition of " << nb_failures
            << " target grid points (source mesh does not contain all target grid points)\n"
            << "Tried first normalizing coordinates with west=" << east - 360.;
        if (second_try) {
            err << "Tried second time normalizing coordinates with west=" << west - eps << "\n";
        }
        err << "Failed target grid points with global index:\n";
        for (size_t n = 0; n < nb_failures; ++n) {
            err << "  - " << std::setw(10) << std::left << failed_index[n] + 1 << " {lon,lat} : " << failed_lonlat[n]
                << "\n";
        }
        throw_Exception(err.str(), Here());
    }
}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
