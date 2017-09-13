/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"

//include <utility>
//include <vector>
//include "eckit/config/Resource.h"
//include "eckit/geometry/Point2.h"
//include "eckit/log/Plural.h"
//include "eckit/log/Timer.h"
//include "eckit/mpi/Comm.h"
#include "eckit/types/FloatCompare.h"
//include "atlas/field/Field.h"
//include "atlas/grid/Grid.h"
//include "atlas/mesh/Elements.h"
//include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
//include "atlas/runtime/Log.h"
//#include "atlas/array.h"
#include "atlas/util/CoordinateEnums.h"


namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {


namespace {
double dot_sign(
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy ) {
  return (Ax - Cx) * (By - Cy)
       - (Bx - Cx) * (Ay - Cy);
}
}


void MatchingMeshPartitioner::getPointCoordinates(const std::vector<idx_t>& poly, std::vector<atlas::PointLonLat>& points, PointLonLat& pointsMin, PointLonLat& pointsMax) const {
    points.clear();
    points.reserve(poly.size());

    auto lonlat = array::make_view< double, 2 >( prePartitionedMesh_.nodes().lonlat() );
    pointsMin = PointLonLat(lonlat(poly[0], LON), lonlat(poly[0], LAT));
    pointsMax = pointsMin;

    for (size_t i = 0; i < poly.size(); ++i) {

        PointLonLat A(lonlat(poly[i], LON), lonlat(poly[i], LAT));
        pointsMin = PointLonLat::componentsMin(pointsMin, A);
        pointsMax = PointLonLat::componentsMax(pointsMax, A);

#if 0
        // if new point is aligned with existing edge (cross product ~= 0) make the edge longer
        if ((points.size() >= 2) && removeAlignedPoints) {
            PointLonLat B = points.back();
            PointLonLat C = points[points.size() - 2];
            if (eckit::types::is_approximately_equal<double>( 0, dot_sign(A[LON], A[LAT], B[LON], B[LAT], C[LON], C[LAT]) )) {
                points.back() = A;
                continue;
            }
        }
#endif

        points.push_back(A);
    }
    ASSERT(points.size() >= 2);
}


}  // partitioner
}  // detail
}  // grid
}  // atlas

