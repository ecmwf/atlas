/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "MatchingMeshPartitionerSphericalPolygon.h"

#include <utility>
#include <vector>
#include "eckit/config/Resource.h"
#include "eckit/log/Plural.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/mpi/Comm.h"
#include "eckit/types/FloatCompare.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"


namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {


PartitionerBuilder<MatchingMeshPartitionerSphericalPolygon> __builder("spherical-polygon");


void MatchingMeshPartitionerSphericalPolygon::partition( const Grid& grid, int node_partition[] ) const {
    eckit::mpi::Comm& comm = eckit::mpi::comm();
    const int mpi_rank = int(comm.rank());

    const Mesh::Polygon& poly = prePartitionedMesh_.polygon(0);

    // Point coordinates (from polygon)
    // - use a bounding box to quickly discard points,
    // - except when that is above/below bounding box but poles should be included
    std::vector< PointLonLat > points;
    PointLonLat pointsMin;
    PointLonLat pointsMax;
    getPointCoordinates(poly, points, pointsMin, pointsMax);
    ASSERT(points.size() >= 2);
    ASSERT(points.size() == poly.size());

    // FIXME: THIS IS A HACK! the coordinates include North/South Pole (first/last partitions only)
    ASSERT( grid.domain().global() );
    bool includes_north_pole = (mpi_rank == 0);
    bool includes_south_pole = (mpi_rank == (int(comm.size()) - 1 ));

    std::vector< PointLonLat > lonlat_tgt_pts;
    lonlat_tgt_pts.reserve(grid.size());
    for( PointXY Pxy : grid.xy() ) {
      lonlat_tgt_pts.push_back( grid.projection().lonlat(Pxy) );
    }

    {
        eckit::ProgressTimer timer("Partitioning target", grid.size(), "point", double(10), atlas::Log::info());
        for (size_t i=0; i<grid.size(); ++i) {
            ++timer;

            node_partition[i] = -1;
            const PointLonLat& P(lonlat_tgt_pts[i]);

            if (pointsMin[LON] <= P[LON] && P[LON] < pointsMax[LON]
             && pointsMin[LAT] <= P[LAT] && P[LAT] < pointsMax[LAT]) {

                if (poly.containsPointInSphericalGeometry(points, P)) {
                    node_partition[i] = mpi_rank;
                }

            } else if ((includes_north_pole && P[LAT] >= pointsMax[LAT])
                    || (includes_south_pole && P[LAT] <  pointsMin[LAT])) {

                node_partition[i] = mpi_rank;

            }
        }
    }


    // Synchronize the partitioning
    comm.allReduceInPlace(node_partition, grid.size(), eckit::mpi::Operation::MAX);


    /// For debugging purposes at the moment. To be made available later, when the Mesh
    /// contains a Polygon for its partition boundary
    if (eckit::Resource<bool>("--output-polygons", false)) {
        poly.outputPythonScript("polygons.py");
    }

    // Sanity check
    const int min = *std::min_element(node_partition, node_partition+grid.size());
    if (min<0) {
        throw eckit::SeriousBug("Could not find partition for input node (meshSource does not contain all points of gridTarget)", Here());
    }
}


}  // partitioner
}  // detail
}  // grid
}  // atlas

