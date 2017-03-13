/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "PrePartitionedBruteForce.h"

#include <vector>
#include "eckit/geometry/Point2.h"
#include "eckit/log/Plural.h"
#include "eckit/log/Timer.h"
#include "eckit/mpi/Comm.h"
#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Log.h"


namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {


namespace {


PartitionerBuilder<PrePartitionedBruteForce> __builder("PrePartitionedBruteForce");


double dot_sign(
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy ) {
  return (Ax - Cx) * (By - Cy)
       - (Bx - Cx) * (Ay - Cy);
}


/// Point-in-triangle test
/// @note "dot product" test, equivalent to half-plane check
/// @see http://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
bool point_in_triangle(
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy,
        const double& Px, const double& Py ) {

    // Compute signs of dot products
    const bool
            b1 = dot_sign(Px, Py, Ax, Ay, Bx, By) < 0,
            b2 = dot_sign(Px, Py, Bx, By, Cx, Cy) < 0,
            b3 = dot_sign(Px, Py, Cx, Cy, Ax, Ay) < 0;

    // Check if point is in triangle
    // - check for b1 && b2 && b3 only works for triangles ordered counter-clockwise
    // - check for b1==b2 && b2==b3 works for both, equivalent to "all points on the same side"
    return (b1 == b2) && (b2 == b3);
}


/// Point-in-quadrilateral test
/// @note limited to convex quadrilaterals, @see https://en.wikipedia.org/wiki/Convex_polygon
bool point_in_quadrilateral (
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy,
        const double& Dx, const double& Dy,
        const double& Px, const double& Py ) {
    const bool
            b1 = dot_sign(Px, Py, Ax, Ay, Bx, By) < 0,
            b2 = dot_sign(Px, Py, Bx, By, Cx, Cy) < 0,
            b3 = dot_sign(Px, Py, Cx, Cy, Dx, Dy) < 0,
            b4 = dot_sign(Px, Py, Dx, Dy, Ax, Ay) < 0;
    return (b1 == b2) && (b2 == b3) && (b3 == b4);
}


enum {LON=0, LAT=1};
typedef eckit::geometry::Point2 point_t;


}  // (anonymous namespace)

void PrePartitionedBruteForce::partition( int node_partition[] ) const {


  eckit::mpi::Comm& comm = eckit::mpi::comm();
  const int mpi_rank = int(comm.rank());


  // Partition bounding box
  ASSERT(prePartitionedMesh_->nodes().size());
  array::ArrayView< double, 2 > lonlat_src( prePartitionedMesh_->nodes().lonlat() );

  point_t bbox_min = point_t(lonlat_src(0, LON), lonlat_src(0, LAT));
  point_t bbox_max = bbox_min;
  for (size_t i = 1; i < prePartitionedMesh_->nodes().size(); ++i) {
      point_t P(lonlat_src(i, LON), lonlat_src(i, LAT));
      bbox_min = point_t::componentsMin(bbox_min, P);
      bbox_max = point_t::componentsMax(bbox_max, P);
  }


  // Partition the target grid nodes
  // - use a polygon bounding box to quickly discard points,
  // - except when that is above/below bounding box but poles should be included
  // std::vector<int> node_partition(grid().npts(), -1);

  for( size_t j=0; j<grid().npts(); ++j )
    node_partition[j] = -1;


  // THIS IS A DIRTY HACK!
  ASSERT( grid().domain().global() );
  bool includes_north_pole = (mpi_rank == 0);
  bool includes_south_pole = (mpi_rank == (int(comm.size()) - 1 ));


  std::vector< PointLonLat > lonlat_tgt_pts;
  lonlat_tgt_pts.reserve(grid().npts());

  for( PointXY Pxy : grid() ) {
    lonlat_tgt_pts.push_back( grid().projection().lonlat(Pxy) );
  }

  {
      eckit::TraceTimer<Atlas> timer("Partitioning target grid...");
      for (size_t i=0; i<grid().npts(); ++i) {

          if (i && (i % 1000 == 0)) {
              double rate = i / timer.elapsed();
              Log::info() << eckit::BigNum(i) << " (at " << rate << " points/s)..." << std::endl;
          }

          const point_t P(lonlat_tgt_pts[i].lon(), lonlat_tgt_pts[i].lat());

          if (bbox_min[LON] <= P[LON] && P[LON] <= bbox_max[LON]
              && bbox_min[LAT] <= P[LAT] && P[LAT] <= bbox_max[LAT]) {

              const mesh::Cells&  elements_src = prePartitionedMesh_->cells();
              const size_t nb_types = elements_src.nb_types();
              bool found = false;

              for (size_t t=0; t<nb_types && !found; ++t) {
                  size_t idx[4];
                  const mesh::Elements& elements = elements_src.elements(t);
                  const mesh::Elements::Connectivity& conn = elements.node_connectivity();

                  const size_t nb_nodes = elements.nb_nodes();
                  ASSERT( (nb_nodes==3 && elements.name() == "Triangle")
                          || (nb_nodes==4 && elements.name() == "Quadrilateral") );

                  for (size_t j=0; j<elements.size() && !found; ++j) {
                      idx[0] = static_cast<size_t>( conn(j,0) );
                      idx[1] = static_cast<size_t>( conn(j,1) );
                      idx[2] = static_cast<size_t>( conn(j,2) );
                      idx[3] = nb_nodes > 3? static_cast<size_t>( conn(j,3) ) : 0;

                      if ((elements.name() == "Triangle" && point_in_triangle(
                               lonlat_src(idx[0], LON), lonlat_src(idx[0], LAT),
                               lonlat_src(idx[1], LON), lonlat_src(idx[1], LAT),
                               lonlat_src(idx[2], LON), lonlat_src(idx[2], LAT),
                               P[LON], P[LAT] ))
                          || (elements.name() == "Quadrilateral" && point_in_quadrilateral(
                                  lonlat_src(idx[0], LON), lonlat_src(idx[0], LAT),
                                  lonlat_src(idx[1], LON), lonlat_src(idx[1], LAT),
                                  lonlat_src(idx[2], LON), lonlat_src(idx[2], LAT),
                                  lonlat_src(idx[3], LON), lonlat_src(idx[3], LAT),
                                  P[LON], P[LAT] )) ) {

                          node_partition[i] = mpi_rank;
                          found = true;

                      }
                  }

              }
          } else if ((includes_north_pole && P[LAT] > bbox_max[LAT])
                  || (includes_south_pole && P[LAT] < bbox_min[LAT])) {

              node_partition[i] = mpi_rank;

          }
      }
  }


  // Synchronize the partitioning and return a grid partitioner
  comm.allReduceInPlace(node_partition, grid().npts(), eckit::mpi::Operation::MAX);
  const int min = *std::min_element(node_partition, node_partition+grid().npts());
  if (min<0) {
      throw eckit::SeriousBug("Could not find partition for input node (partitionedMesh does not contain all points of gridToDistribute)", Here());
  }

}


}  // partitioner
}  // detail
}  // grid
}  // atlas

