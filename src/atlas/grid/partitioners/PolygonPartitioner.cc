/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "PolygonPartitioner.h"

#include <set>
#include <utility>
#include <vector>
#include "eckit/geometry/Point2.h"
#include "eckit/log/Plural.h"
#include "eckit/log/Timer.h"
#include "eckit/mpi/Comm.h"
#include "eckit/types/FloatCompare.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Structured.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/LibAtlas.h"
#include "atlas/runtime/Log.h"


namespace atlas {
namespace grid {
namespace partitioners {


namespace {


PartitionerBuilder<PolygonPartitioner> __builder("polygon");


double dot_sign(
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy ) {
  return (Ax - Cx) * (By - Cy)
       - (Bx - Cx) * (Ay - Cy);
}


enum {LON=0, LAT=1};
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


struct Edge : std::pair< size_t, size_t > {
    Edge(size_t _A, size_t _B) : pair(_A, _B) {}
    Edge reverse() const { return Edge(second, first); }
    struct LessThan {
        bool operator()(const Edge& e1, const Edge& e2) const {
            // order ascending by 'first'
            return (e1.first  < e2.first? true  :
                    e1.first  > e2.first? false :
                    e1.second < e2.second );
        }
    };
};


typedef std::set< Edge, Edge::LessThan > edge_set_t;


struct Poly : std::vector< size_t > {
    Poly() {}
    Poly(edge_set_t& edges) {
        clear();
        reserve(edges.size() + 1);

        push_back(edges.begin()->first);
        for (edge_set_t::iterator e = edges.begin(); e != edges.end() && e->first == back();
             e = edges.lower_bound(Edge(back(), 0))) {
            push_back(e->second);
            edges.erase(*e);
        }

        ASSERT(front() == back());
    }
    Poly& operator+=(const Poly& other) {
        ASSERT(other.size() > 1);
        if (empty()) {
            return operator=(other);
        }

        difference_type N = difference_type(other.size()) - 1;
        vector< size_t > cycle;
        cycle.reserve(2 * size_t(N));
        cycle.insert(cycle.end(), other.begin(), other.end() - 1);
        cycle.insert(cycle.end(), other.begin(), other.end() - 1);

        for (const_iterator c = cycle.begin(); c != cycle.begin() + N; ++c) {
            iterator here = find(begin(), end(), *c);
            if (here != end()) {
                insert(here, c, c + N);
                return *this;
            }
        }

        throw eckit::AssertionFailed("Could not merge polygons, they are not connected", Here());
    }
};


}  // (anonymous namespace)


GridDistribution* PolygonPartitioner::distributionFromPrePartitionedMesh() const {
    eckit::mpi::Comm& comm = eckit::mpi::comm();
    const int mpi_rank = int(comm.rank());


    // Polygon indices
    // - extract partition boundary edges by always attempting first to
    //   remove a reversed edge from a neighbouring element (if it exists)
    // - polygon can have multiple cycles, but must be a connected graph
    Poly poly;

    edge_set_t edges;
    for (size_t t = 0; t < prePartitionedMesh_->cells().nb_types(); ++t) {
        const mesh::Elements& elements = prePartitionedMesh_->cells().elements(t);

        const mesh::Elements::Connectivity& conn = elements.node_connectivity();
        array::ArrayView< int, 1 > patch( elements.view< int, 1 >(elements.field("patch")) );

        const size_t nb_nodes = elements.nb_nodes();
        ASSERT( (nb_nodes==3 && elements.name() == "Triangle")
             || (nb_nodes==4 && elements.name() == "Quadrilateral") );

        for (size_t j = 0; j < elements.size(); ++j) {
            for (size_t k = 0; k < nb_nodes && (patch(j) == 0); ++k) {
                Edge edge(size_t(conn(j, k)), size_t(conn(j, (k+1) % nb_nodes)));
                ASSERT(edge.first != edge.second);
                if (!edges.erase(edge.reverse())) {
                    edges.insert(edge);
                }
            }
        }
    }

    while (!edges.empty()) {
        poly += Poly(edges);
    }
    ASSERT(poly.size());


    // Polygon point coordinates: simplify by checking aligned edges
    // Note: indices ('poly') don't necessarily match coordinates ('polygon')
    // Note: the coordinates include North/South Pole (first/last partitions only)
    std::vector< point_t > polygon;
    polygon.reserve(polygon.size());

    array::ArrayView< double, 2 > lonlat_src( prePartitionedMesh_->nodes().lonlat() );
    point_t bbox_min = point_t(lonlat_src(poly[0], LON), lonlat_src(poly[0], LAT));
    point_t bbox_max = bbox_min;

    for (size_t i = 0; i < poly.size(); ++i) {

        point_t A(lonlat_src(poly[i], LON), lonlat_src(poly[i], LAT));
        bbox_min = point_t::componentsMin(bbox_min, A);
        bbox_max = point_t::componentsMax(bbox_max, A);

        // if new point is aligned with existing edge (cross product ~= 0) make the edge longer
        if ((polygon.size() >= 2)) {
            point_t B = polygon.back();
            point_t C = polygon[polygon.size() - 2];
            if (eckit::types::is_approximately_equal<double>( 0, dot_sign(A[LON], A[LAT], B[LON], B[LAT], C[LON], C[LAT]) )) {
                polygon.back() = A;
                continue;
            }
        }

        polygon.push_back(A);
    }
    ASSERT(polygon.size() >= 2);


    // Partition the target grid nodes
    // - use a polygon bounding box to quickly discard points,
    // - except when that is above/below bounding box but poles should be included
    std::vector<int> node_partition(grid().npts(), -1);

    bool includes_north_pole = prePartitionedDomain_.includesPoleNorth() && mpi_rank == 0;
    bool includes_south_pole = prePartitionedDomain_.includesPoleSouth() && mpi_rank == (int(comm.size()) - 1 );

    std::vector< grid::Grid::Point > lonlat_tgt_pts;
    grid().lonlat(lonlat_tgt_pts);
    {
        eckit::TraceTimer<LibAtlas> timer("Partitioning target grid...");
        for (size_t i=0; i<grid().npts(); ++i) {

            if (i && (i % 1000 == 0)) {
                double rate = i / timer.elapsed();
                Log::info() << eckit::BigNum(i) << " (at " << rate << " points/s)..." << std::endl;
            }

            point_t P(lonlat_tgt_pts[i].lon(), lonlat_tgt_pts[i].lat());

            if (bbox_min[LON] <= P[LON] && P[LON] <= bbox_max[LON]
             && bbox_min[LAT] <= P[LAT] && P[LAT] <= bbox_max[LAT]) {

                if (point_in_poly(polygon, P)) {
                    node_partition[i] = mpi_rank;
                }

            } else if ((includes_north_pole && P[LAT] > bbox_max[LAT])
                    || (includes_south_pole && P[LAT] < bbox_min[LAT])) {

                node_partition[i] = mpi_rank;

            }
        }
    }


    // Synchronize the partitioning and return a grid partitioner
    comm.allReduceInPlace(node_partition.data(), node_partition.size(), eckit::mpi::Operation::MAX);
    const int min = *std::min_element(node_partition.begin(), node_partition.end());
    if (min<0) {
        throw eckit::SeriousBug("Could not find partition for input node (meshSource does not contain all points of gridTarget)", Here());
    }

    GridDistribution* distribution_tgt = new grid::GridDistribution(node_partition.size(), node_partition.data());
    comm.barrier();
#if 0
    std::vector<double> x,y;
    x.reserve(node_partition.size());
    y.reserve(node_partition.size());
    for (size_t i=0; i<node_partition.size(); ++i) {
        if (node_partition[i] == mpi_rank) {
            x.push_back(lonlat_tgt_pts[i].lon());
            y.push_back(lonlat_tgt_pts[i].lat());
        }
    }
    size_t count = x.size();
    size_t count_all = x.size();
    comm.allReduceInPlace(count_all, eckit::mpi::sum());

    for (int r = 0; r < comm.size(); ++r) {
        if (mpi_rank == r) {
            std::ofstream f("partitions_poly.py", mpi_rank == 0? std::ios::trunc : std::ios::app);

            if (mpi_rank == 0) {
                f << "\n" "import matplotlib.pyplot as plt"
                     "\n" "from matplotlib.path import Path"
                     "\n" "import matplotlib.patches as patches"
                     "\n" ""
                     "\n" "from itertools import cycle"
                     "\n" "cycol = cycle('bgrcmy').next"
                     "\n" ""
                     "\n" "fig = plt.figure()"
                     "\n" "ax = fig.add_subplot(111)"
                     "\n" "";
            }
            f << "\n" "verts_" << r << " = [";
            for (size_t i = 0; i < poly.size(); ++i) { f << "\n  (" << poly[i][LON] << ", " << poly[i][LAT] << "), "; }
            f << "\n]"
                 "\n" ""
                 "\n" "codes_" << r << " = [Path.MOVETO]"
                 "\n" "codes_" << r << ".extend([Path.LINETO] * " << (poly.size()-2) << ")"
                 "\n" "codes_" << r << ".extend([Path.CLOSEPOLY])"
                 "\n" ""
                 "\n" "count_" << r << " = " << count <<
                 "\n" "count_all_" << r << " = " << count_all <<
                 "\n" ""
                 "\n" "x_" << r << " = ["; for (std::vector<double>::const_iterator ix=x.begin(); ix!=x.end(); ++ix) { f << *ix << ", "; } f << "]"
                 "\n" "y_" << r << " = ["; for (std::vector<double>::const_iterator iy=y.begin(); iy!=y.end(); ++iy) { f << *iy << ", "; } f << "]"
                 "\n"
                 "\n" "c = cycol()"
                 "\n" "ax.add_patch(patches.PathPatch(Path(verts_" << r << ", codes_" << r << "), facecolor=c, color=c, alpha=0.3, lw=1))"
                 "\n" "ax.scatter(x_" << r << ", y_" << r << ", color=c, marker='o')"
                 "\n" "";
            if (mpi_rank == int(comm.size()) - 1) {
                f << "\n" "ax.set_xlim(  0-5, 360+5)"
                     "\n" "ax.set_ylim(-90-5,  90+5)"
                     "\n" "plt.show()";
            }
        }
        comm.barrier();
    }
//    exit(0);
#endif
    return distribution_tgt;
}


}  // partitioners
}  // grid
}  // atlas

