/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerLonLatPolygon.h"

#include <vector>
#include "eckit/config/Resource.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/mpi/Comm.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Log.h"


namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {


namespace {
PartitionerBuilder<MatchingMeshPartitionerLonLatPolygon> __builder("lonlat-polygon");
}


void MatchingMeshPartitionerLonLatPolygon::partition( const Grid& grid, int partitioning[] ) const {
    eckit::mpi::Comm& comm = eckit::mpi::comm();
    const int mpi_rank = int(comm.rank());

    ASSERT( grid.domain().global() );
    const Mesh::Polygon& poly = prePartitionedMesh_.polygon(0);
    poly.setCoordinates();

    {
        eckit::ProgressTimer timer("Partitioning target", grid.size(), "point", double(10), atlas::Log::info());
        size_t i = 0;

        for (PointXY Pxy : grid.xy()) {
            ++timer;
            const PointLonLat P = grid.projection().lonlat(Pxy);
            partitioning[i++] = poly.containsPointInLonLatGeometry(P) ? mpi_rank : -1;
        }

        // Synchronize partitioning, do a sanity check
        comm.allReduceInPlace(partitioning, grid.size(), eckit::mpi::Operation::MAX);
        const int min = *std::min_element(partitioning, partitioning+grid.size());
        if (min<0) {
            throw eckit::SeriousBug("Could not find partition for target node (source mesh does not contain all target grid points)", Here());
        }
    }

    /// For debugging purposes
    if (eckit::Resource<bool>("--output-polygons", false)) {
        poly.outputPythonScript("polygons.py");
    }
}


}  // partitioner
}  // detail
}  // grid
}  // atlas

