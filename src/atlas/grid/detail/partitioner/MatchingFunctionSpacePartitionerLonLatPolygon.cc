/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/partitioner/MatchingFunctionSpacePartitionerLonLatPolygon.h"

#include <vector>

#include "eckit/config/Resource.h"
#include "eckit/log/ProgressTimer.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/PolygonXY.h"

#include "atlas/parallel/omp/fill.h"
#include "atlas/parallel/omp/omp.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

namespace {
// We did not yet implement self-registration in MatchedPartitionerFactory
// PartitionerBuilder<MatchingFunctionSpacePartitionerLonLatPolygon> __builder( "lonlat-polygon" );
}

void MatchingFunctionSpacePartitionerLonLatPolygon::partition(const Grid& grid, int part[]) const {
    ATLAS_TRACE("MatchingFunctionSpacePartitionerLonLatPolygon");
    //atlas::vector<int> part( grid.size() );

    const auto& comm   = mpi::comm(partitioned_.mpi_comm());
    const int mpi_rank = int(comm.rank());
    const int mpi_size = int(comm.size());

    if (mpi_size == 1) {
        // shortcut
        omp::fill(part, part + grid.size(), 0);
    }
    else {
        const auto& p = partitioned_.polygon();

        util::PolygonXY poly{p};
        {
            ATLAS_TRACE("point-in-polygon check for entire grid (" + std::to_string(grid.size()) + " points)");
            size_t num_threads = atlas_omp_get_max_threads();
            size_t chunk_size  = grid.size() / (1000 * num_threads);
            size_t chunks      = num_threads == 1 ? 1 : std::max(size_t(1), size_t(grid.size()) / chunk_size);
            atlas_omp_pragma(omp parallel for schedule(dynamic,1))
            for( size_t chunk=0; chunk < chunks; ++chunk) {
                const size_t begin = chunk * size_t(grid.size()) / chunks;
                const size_t end   = (chunk + 1) * size_t(grid.size()) / chunks;
                auto it            = grid.xy().begin() + chunk * grid.size() / chunks;
                for (size_t n = begin; n < end; ++n) {
                    if (poly.contains(*it)) {
                        part[n] = mpi_rank;
                    }
                    else {
                        part[n] = -1;
                    }
                    ++it;
                }
            }
        }
        ATLAS_TRACE_MPI(ALLREDUCE) { comm.allReduceInPlace(part, grid.size(), eckit::mpi::max()); }
    }
}


#if 0
    const auto& comm   = mpi::comm(partitioned_.mpi_comm());
    const int mpi_rank = int(comm.rank());
    const int mpi_size = int(comm.size());

    ATLAS_TRACE( "MatchingFunctionSpacePartitionerLonLatPolygon::partition" );

    ATLAS_ASSERT( grid.domain().global() );

    Log::debug() << "MatchingFunctionSpacePartitionerLonLatPolygon::partition" << std::endl;

    // FIXME: THIS IS A HACK! the coordinates include North/South Pole (first/last
    // partitions only)
    bool includesNorthPole = ( mpi_rank == 0 );
    bool includesSouthPole = ( mpi_rank == mpi_size - 1 );

    const util::PolygonXY poly{prePartitionedFunctionSpace_.polygon( 0 )};

    {
        eckit::ProgressTimer timer( "Partitioning", grid.size(), "point", double( 10 ), atlas::Log::trace() );
        size_t i = 0;

        for ( const PointLonLat& P : grid.lonlat() ) {
            ++timer;
            const bool atThePole = ( includesNorthPole && P[LAT] >= poly.coordinatesMax()[LAT] ) ||
                                   ( includesSouthPole && P[LAT] < poly.coordinatesMin()[LAT] );

            partitioning[i++] = atThePole || poly.contains( P ) ? mpi_rank : -1;
        }
    }

    // Synchronize partitioning, do a sanity check
    comm.allReduceInPlace( partitioning, grid.size(), eckit::mpi::Operation::MAX );
    const int min = *std::min_element( partitioning, partitioning + grid.size() );
    if ( min < 0 ) {
        throw_Exception(
            "Could not find partition for target node (source "
            "FunctionSpace does not contain all target grid points)",
            Here() );
    }
#endif

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
