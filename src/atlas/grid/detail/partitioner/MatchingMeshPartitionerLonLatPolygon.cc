/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "MatchingMeshPartitionerLonLatPolygon.h"

#include <utility>
#include <vector>
#include "eckit/config/Resource.h"
#include "eckit/geometry/Point2.h"
#include "eckit/log/Plural.h"
#include "eckit/log/Timer.h"
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


PartitionerBuilder<MatchingMeshPartitionerLonLatPolygon> __builder("lonlat-polygon");


double dot_sign(
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy ) {
  return (Ax - Cx) * (By - Cy)
       - (Bx - Cx) * (Ay - Cy);
}


typedef eckit::geometry::Point2 point_t;


/*
 * Tests if a given point is left|on|right of an infinite line.
 * @input P point to test
 * @input A, B points on infinite line
 * @return >0/=0/<0 for P left|on|right of the infinite line
 */
double point_on_which_side(const point_t& P, const point_t& A, const point_t& B) {
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
bool point_in_poly(const std::vector<point_t>& poly, const point_t& P) {
    ASSERT(poly.size());

    // winding number
    int wn = 0;

    // loop on polygon edges
    for (size_t i = 1; i < poly.size(); ++i) {
        const point_t& A = poly[i-1];
        const point_t& B = poly[ i ];

        // intersect either:
        // - "up" on upward crossing & P left of edge, or
        // - "down" on downward crossing & P right of edge
        if (A[LAT] <= P[LAT] && P[LAT] < B[LAT]) {
            if (point_on_which_side(P, A, B) > 0) {
                ++wn;
            }
        } else if (B[LAT] <= P[LAT] && P[LAT] < A[LAT]) {
            if (point_on_which_side(P, A, B) < 0) {
                --wn;
            }
        }
    }

    // wn == 0 only when P is outside
    return wn != 0;
}

}  // (anonymous namespace)


void MatchingMeshPartitionerLonLatPolygon::partition( const Grid& grid, int node_partition[] ) const {
    eckit::mpi::Comm& comm = eckit::mpi::comm();
    const int mpi_rank = int(comm.rank());
    const int mpi_size = int(comm.size());


    const Mesh::Polygon& _poly = prePartitionedMesh_.polygon(0);

    // Polygon point coordinates: simplify by checking aligned edges
    // Note: indices ('poly') don't necessarily match coordinates ('polygon')
    // Note: the coordinates include North/South Pole (first/last partitions only)
    std::vector< point_t > polygon;
    getPointCoordinates(grid, _poly, polygon);
    ASSERT(polygon.size() >= 2);


    // Partition the target grid nodes
    // - use a polygon bounding box to quickly discard points,
    // - except when that is above/below bounding box but poles should be included
    for( size_t j=0; j<grid.size(); ++j ) {
      node_partition[j] = -1;
    }

    // THIS IS A DIRTY HACK!
    ASSERT( grid.domain().global() );
    bool includes_north_pole = (mpi_rank == 0);
    bool includes_south_pole = (mpi_rank == (int(comm.size()) - 1 ));

    std::vector< PointLonLat > lonlat_tgt_pts;
    lonlat_tgt_pts.reserve(grid.size());

    for( PointXY Pxy : grid.xy() ) {
      lonlat_tgt_pts.push_back( grid.projection().lonlat(Pxy) );
    }

    {
        std::stringstream msg; msg << "Partitioning " << eckit::BigNum(grid.size())
          << " target grid points... ";
        Log::debug<Atlas>() << msg.str() << std::endl;
        eckit::TraceTimer<Atlas> timer(msg.str()+"done");
        for (size_t i=0; i<grid.size(); ++i) {

            if (i && (i % 1000 == 0)) {
                double rate = i / timer.elapsed();
                Log::debug<Atlas>() << "    " << eckit::BigNum(i) << " points completed (at " << rate << " points/s)" << std::endl;
            }

            point_t P(lonlat_tgt_pts[i].lon(), lonlat_tgt_pts[i].lat());

            if (bbox_min[LON] <= P[LON] && P[LON] < bbox_max[LON]
             && bbox_min[LAT] <= P[LAT] && P[LAT] < bbox_max[LAT]) {

                if (point_in_poly(polygon, P)) {
                    node_partition[i] = mpi_rank;
                }

            } else if ((includes_north_pole && P[LAT] >= bbox_max[LAT])
                    || (includes_south_pole && P[LAT] <  bbox_min[LAT])) {

                node_partition[i] = mpi_rank;

            }
        }
    }


    // Synchronize the partitioning and return a grid partitioner
    comm.allReduceInPlace(node_partition, grid.size(), eckit::mpi::Operation::MAX);


    /// For debugging purposes at the moment. To be made available later, when the Mesh
    /// contains a Polygon for its partition boundary
    if( eckit::Resource<bool>("--polygons",false) ) {

        std::vector<double> x,y, xlost,ylost;
        xlost.reserve(grid.size());
        ylost.reserve(grid.size());
        x.reserve(grid.size());
        y.reserve(grid.size());
        for (size_t i=0; i<grid.size(); ++i) {
            if (node_partition[i] == mpi_rank) {
                x.push_back(lonlat_tgt_pts[i].lon());
                y.push_back(lonlat_tgt_pts[i].lat());
            } else if (node_partition[i] == -1) {
                xlost.push_back(lonlat_tgt_pts[i].lon());
                ylost.push_back(lonlat_tgt_pts[i].lat());
            }

        }
        size_t count = x.size();
        size_t count_all = x.size();
        comm.allReduceInPlace(count_all, eckit::mpi::sum());

        for (int r = 0; r < int(comm.size()); ++r) {
            if (mpi_rank == r) {
                std::ofstream f("partitions_poly.py", mpi_rank == 0? std::ios::trunc : std::ios::app);

                if (mpi_rank == 0) {
                    f << "\n" "import matplotlib.pyplot as plt"
                         "\n" "from matplotlib.path import Path"
                         "\n" "import matplotlib.patches as patches"
                         "\n" ""
                         "\n" "from itertools import cycle"
                         "\n" "import matplotlib.cm as cm"
                         "\n" "import numpy as np"
                         "\n" "cycol = cycle([cm.Paired(i) for i in np.linspace(0,1,12,endpoint=True)]).next"
                         "\n" ""
                         "\n" "fig = plt.figure()"
                         "\n" "ax = fig.add_subplot(111,aspect='equal')"
                         "\n" "";
                    f << "\n"
                         "\n" "xlost = ["; for (const double& ix: xlost) { f << ix << ", "; } f << "]"
                         "\n" "ylost = ["; for (const double& iy: ylost) { f << iy << ", "; } f << "]"
                         "\n"
                         "\n" "ax.scatter(xlost, ylost, color='k', marker='o')"
                         "\n" "";
                }
                f << "\n" "verts_" << r << " = [";
                for (size_t i = 0; i < polygon.size(); ++i) { f << "\n  (" << polygon[i][LON] << ", " << polygon[i][LAT] << "), "; }
                f << "\n]"
                     "\n" ""
                     "\n" "codes_" << r << " = [Path.MOVETO]"
                     "\n" "codes_" << r << ".extend([Path.LINETO] * " << (polygon.size()-2) << ")"
                     "\n" "codes_" << r << ".extend([Path.CLOSEPOLY])"
                     "\n" ""
                     "\n" "count_" << r << " = " << count <<
                     "\n" "count_all_" << r << " = " << count_all <<
                     "\n" ""
                     "\n" "x_" << r << " = ["; for (const double& ix: x) { f << ix << ", "; } f << "]"
                     "\n" "y_" << r << " = ["; for (const double& iy: y) { f << iy << ", "; } f << "]"
                     "\n"
                     "\n" "c = cycol()"
                     "\n" "ax.add_patch(patches.PathPatch(Path(verts_" << r << ", codes_" << r << "), facecolor=c, color=c, alpha=0.3, lw=1))"
                     "\n" "ax.scatter(x_" << r << ", y_" << r << ", color=c, marker='o')"
                     "\n" "";
                if (mpi_rank == int(comm.size()) - 1) {
                    f << "\n" "ax.set_xlim(  0-5, 360+5)"
                         "\n" "ax.set_ylim(-90-5,  90+5)"
                         "\n" "ax.set_xticks([0,45,90,135,180,225,270,315,360])"
                         "\n" "ax.set_yticks([-90,-45,0,45,90])"
                         "\n" "plt.grid()"
                         "\n" "plt.show()";
                }
            }
            comm.barrier();
        }

    }


    // Sanity check
    {
      const int min = *std::min_element(node_partition, node_partition+grid.size());
      if (min<0) {
          throw eckit::SeriousBug("Could not find partition for input node (meshSource does not contain all points of gridTarget)", Here());
      }
    }

}


void MatchingMeshPartitionerLonLatPolygon::getPointCoordinates(const Grid&, const Mesh::Polygon& poly, std::vector<PointLonLat>& points) const {
    points.clear();
    points.reserve(poly.size());

    auto lonlat = array::make_view< double, 2 >( prePartitionedMesh_.nodes().lonlat() );
    PointLonLat bbox_min = PointLonLat(lonlat(poly[0], LON), lonlat(poly[0], LAT));
    PointLonLat bbox_max = bbox_min;

    for (size_t i = 0; i < poly.size(); ++i) {

        PointLonLat A(lonlat(poly[i], LON), lonlat(poly[i], LAT));
        bbox_min = PointLonLat::componentsMin(bbox_min, A);
        bbox_max = PointLonLat::componentsMax(bbox_max, A);

        // if new point is aligned with existing edge (cross product ~= 0) make the edge longer
        if ((points.size() >= 2)) {
            PointLonLat B = points.back();
            PointLonLat C = points[points.size() - 2];
            if (eckit::types::is_approximately_equal<double>( 0, dot_sign(A[LON], A[LAT], B[LON], B[LAT], C[LON], C[LAT]) )) {
                points.back() = A;
                continue;
            }
        }

        points.push_back(A);
    }
    ASSERT(points.size() >= 2);
}


}  // partitioner
}  // detail
}  // grid
}  // atlas

