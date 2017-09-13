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


namespace {


PartitionerBuilder<MatchingMeshPartitionerSphericalPolygon> __builder("spherical-polygon");


double dot_sign(
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy ) {
  return (Ax - Cx) * (By - Cy)
       - (Bx - Cx) * (Ay - Cy);
}


/*
 * Tests if a given point is left|on|right of an infinite line.
 * @input P point to test
 * @input A, B points on infinite line
 * @return >0/=0/<0 for P left|on|right of the infinite line
 */
double point_on_which_side(const PointLonLat& P, const PointLonLat& A, const PointLonLat& B) {
    return dot_sign( P[LON], P[LAT],
                     A[LON], A[LAT],
                     B[LON], B[LAT] );
}


/*
 * Point inclusion based on winding number for a point in a polygon
 * @note reference <a href="http://geomalgorithms.com/a03-_inclusion.html">Inclusion of a Point in a Polygon</a>
 * @param[in] poly vertex points of a polygon (closed, where poly.front() == poly.back())
 * @param[in] P given point
 * @return if point is in polygon
 */
bool point_in_poly(const std::vector<PointLonLat>& poly, const PointLonLat& P) {
    ASSERT(poly.size());

    // winding number
    int wn = 0;

    // loop on polygon edges
    for (size_t i = 1; i < poly.size(); ++i) {
        const PointLonLat& A = poly[i-1];
        const PointLonLat& B = poly[ i ];

        // intersect either:
        // - "up" on upward crossing & P left of edge, or
        // - "down" on downward crossing & P right of edge
        if (A[LON] <= P[LON] && P[LON] < B[LON]) {
            if (point_on_which_side(P, A, B) > 0) {
                ++wn;
            }
        } else if (B[LON] <= P[LON] && P[LON] < A[LON]) {
            if (point_on_which_side(P, A, B) < 0) {
                --wn;
            }
        }
    }

    // wn == 0 only when P is outside
    return wn != 0;
}

}  // (anonymous namespace)


void MatchingMeshPartitionerSphericalPolygon::partition( const Grid& grid, int node_partition[] ) const {
    eckit::mpi::Comm& comm = eckit::mpi::comm();
    const int mpi_rank = int(comm.rank());

    const Mesh::Polygon& _poly = prePartitionedMesh_.polygon(0);

    // Point coordinates (from polygon)
    // - use a bounding box to quickly discard points,
    // - except when that is above/below bounding box but poles should be included
    std::vector< PointLonLat > polygon;
    PointLonLat pointsMin;
    PointLonLat pointsMax;
    getPointCoordinates(_poly, polygon, pointsMin, pointsMax, false);
    ASSERT(polygon.size() >= 2);
    ASSERT(polygon.size() == _poly.size());

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

                if (point_in_poly(polygon, P)) {
                    node_partition[i] = mpi_rank;
                }

            } else if ((includes_north_pole && P[LAT] >= pointsMax[LAT])
                    || (includes_south_pole && P[LAT] <  pointsMin[LAT])) {

                node_partition[i] = mpi_rank;

            }
        }
    }


    // Synchronize the partitioning and return a grid partitioner
    comm.allReduceInPlace(node_partition, grid.size(), eckit::mpi::Operation::MAX);


    /// For debugging purposes at the moment. To be made available later, when the Mesh
    /// contains a Polygon for its partition boundary
    if (eckit::Resource<bool>("--polygons", false)) {
        _poly.outputPythonScript("partitions_poly.py");
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

