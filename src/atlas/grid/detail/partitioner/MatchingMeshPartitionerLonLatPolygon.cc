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
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/PolygonXY.h"


#include "atlas/util/KDTree.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

namespace {
PartitionerBuilder<MatchingMeshPartitionerLonLatPolygon> __builder("lonlat-polygon");
}

MatchingMeshPartitionerLonLatPolygon::MatchingMeshPartitionerLonLatPolygon(const Mesh& mesh, const eckit::Parametrisation& config):
    MatchingMeshPartitioner(mesh, config) {
        config.get("fallback_nearest", fallback_nearest_);
    }


void MatchingMeshPartitionerLonLatPolygon::partition(const Grid& grid, int partitioning[]) const {
    const auto& comm   = mpi::comm(prePartitionedMesh_.mpi_comm());
    const int mpi_rank = int(comm.rank());
    const int mpi_size = int(comm.size());

    ATLAS_TRACE("MatchingMeshPartitionerLonLatPolygon::partition");

    ATLAS_ASSERT(grid.domain().global());

    Log::debug() << "MatchingMeshPartitionerLonLatPolygon::partition" << std::endl;

    const util::PolygonXY poly{prePartitionedMesh_.polygon(0)};

    double west = poly.coordinatesMin().x();
    double east = poly.coordinatesMax().x();
    ATLAS_TRACE_MPI(ALLREDUCE) {
      comm.allReduceInPlace(west, eckit::mpi::Operation::MIN);
      comm.allReduceInPlace(east, eckit::mpi::Operation::MAX);
    }

    Projection projection = prePartitionedMesh_.projection();
    omp::fill(partitioning, partitioning + grid.size(), -1);

    auto compute = [&](double west) {
#if !defined(__NVCOMPILER)
        atlas_omp_parallel
#elif (__NVCOMPILER_MAJOR__ >= 21) && (__NVCOMPILER_MINOR__ > 9 )
        // Internal compiler error with nvhpc 21.9:
        //
        //    NVC++-S-0000-Internal compiler error. BAD sptr in var_refsym       0  (MatchingMeshPartitionerLonLatPolygon.cc: 64) (=following line)
        //    NVC++/x86-64 Linux 21.9-0: compilation completed with severe errors
        atlas_omp_parallel
#endif
        {
            const idx_t num_threads = atlas_omp_get_num_threads();
            const idx_t thread_num  = atlas_omp_get_thread_num();
            const idx_t begin = static_cast<idx_t>(thread_num * size_t(grid.size()) / num_threads);
            const idx_t end =
                static_cast<idx_t>((thread_num + 1) * size_t(grid.size()) / num_threads);
            size_t i = begin;
            auto it = grid.lonlat().begin() + i;
            for(; i < end; ++i, ++it) {
                PointLonLat P = *it;
                if (partitioning[i] < 0) {
                    projection.lonlat2xy(P);
                    P.normalise(west);
                    partitioning[i] = poly.contains(P) ? mpi_rank : -1;
                }
            }
        }
        // Synchronize partitioning
        ATLAS_TRACE_MPI(ALLREDUCE) {
          comm.allReduceInPlace(partitioning, grid.size(), eckit::mpi::Operation::MAX);
        }

        std::vector<int> thread_min(atlas_omp_get_max_threads(),std::numeric_limits<int>::max());
#if !defined(__NVCOMPILER)
        atlas_omp_parallel
#elif (__NVCOMPILER_MAJOR__ >= 21) && (__NVCOMPILER_MINOR__ > 9 )
        // Internal compiler error with nvhpc 21.9:
        //
        //    NVC++-S-0000-Internal compiler error. BAD sptr in var_refsym       0  (MatchingMeshPartitionerLonLatPolygon.cc: 64) (=following line)
        //    NVC++/x86-64 Linux 21.9-0: compilation completed with severe errors
        atlas_omp_parallel
#endif
       {
            const idx_t num_threads = atlas_omp_get_num_threads();
            const idx_t thread_num  = atlas_omp_get_thread_num();
            const idx_t begin = static_cast<idx_t>(thread_num * size_t(grid.size()) / num_threads);
            const idx_t end =
                static_cast<idx_t>((thread_num + 1) * size_t(grid.size()) / num_threads);
            thread_min[thread_num] = *std::min_element(partitioning+begin,partitioning+end);
        }
        return *std::min_element(thread_min.begin(), thread_min.end());
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

        if (fallback_nearest_) {
            util::IndexKDTree3D kdtree;
            kdtree.reserve(grid.size());
            size_t j=0;
            for (auto& p: grid.lonlat()) {
                if (partitioning[j] >= 0) {
                    kdtree.insert(p,partitioning[j]);
                }
                ++j;
            }
            kdtree.build();
            for (size_t n = 0; n < nb_failures; ++n) {
                auto closest = kdtree.closestPoint(failed_lonlat[n]);
                partitioning[failed_index[n]] = closest.payload();
            }
        }
        else {
            std::stringstream err;
            err.precision(20);
            err << "Could not find partition of " << nb_failures
                << " target grid points (source mesh does not contain all target grid points)\n"
                << "Tried first normalizing coordinates with west=" << east - 360. << "\n";
            if (second_try) {
                err << "Tried second time normalizing coordinates with west=" << west - eps << "\n";
            }

            err << "Failed target grid points with global index:\n";
            for (size_t n = 0; n < nb_failures; ++n) {
                err << "  - " << std::setw(10) << std::left << failed_index[n] + 1 << " {lon,lat} : " << failed_lonlat[n]
                    << "\n";
            }
            // prePartitionedMesh_.polygon(0).outputPythonScript("partitions.py");
            throw_Exception(err.str(), Here());
        }
    }
}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
