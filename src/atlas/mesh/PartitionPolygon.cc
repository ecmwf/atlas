/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/PartitionPolygon.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
//#include "eckit/config/Resource.h"
//#include "eckit/mpi/Comm.h"
//#include "eckit/types/FloatCompare.h"
//#include "atlas/grid/Grid.h"
//#include "atlas/mesh/Elements.h"
//#include "atlas/mesh/Mesh.h"
//#include "atlas/mesh/Nodes.h"
//#include "atlas/runtime/Log.h"

namespace atlas {
namespace mesh {


//------------------------------------------------------------------------------------------------------


namespace {


double dot_sign(
        const double& Ax, const double& Ay,
        const double& Bx, const double& By,
        const double& Cx, const double& Cy ) {
  return (Ax - Cx) * (By - Cy)
       - (Ay - Cy) * (Bx - Cx);
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


}  // (anonymous)


//------------------------------------------------------------------------------------------------------


PartitionPolygon::PartitionPolygon(const detail::MeshImpl& mesh, size_t halo) :
  mesh_(mesh),
  halo_(halo) {
  const eckit::mpi::Comm& comm = parallel::mpi::comm();
  const int mpi_size = int(comm.size());

  // extract partition boundary edges by always attempting first to`
  // remove a reversed edge from a neighbouring element, if any
  edge_set_t edges;
  for (size_t t = 0; t < mesh.cells().nb_types(); ++t) {
      const Elements& elements = mesh.cells().elements(t);

      const BlockConnectivity& conn = elements.node_connectivity();
      auto field_patch = elements.view< int, 1 >(elements.field("patch"));
      auto field_halo  = elements.view< int, 1 >(elements.field("halo"));

      const size_t nb_nodes = elements.nb_nodes();

      for (size_t j = 0; j < elements.size(); ++j) {
          if (field_patch(j) == 0 && field_halo(j) <= halo ) {
              for (size_t k = 0; k < nb_nodes; ++k) {
                  edge_t edge(size_t(conn(j, k)), size_t(conn(j, (k+1) % nb_nodes)));
                  if (!edges.erase(edge.reverse())) {
                      edges.insert(edge);
                  }
              }
          }
      }
  }

  ASSERT(operator+=(detail::Polygon(edges)));
}


size_t PartitionPolygon::footprint() const {
    size_t size = sizeof(*this);
    size += capacity()*sizeof(idx_t);
    return size;
}


void PartitionPolygon::print(std::ostream& out) const {
    out << "polygon:{"
        <<  "halo:" << halo_
        << ",size:" << size()
        << ",nodes:" << static_cast<detail::Polygon>(*this)
        << "}";
}


void PartitionPolygon::outputPythonScript(const eckit::PathName& filepath) const {
  const eckit::mpi::Comm& comm = atlas::parallel::mpi::comm();
  int mpi_rank = comm.rank();
  int mpi_size = comm.size();

  auto xy = array::make_view<double,2>( mesh_.nodes().xy() );

  double xmin =  std::numeric_limits<double>::max();
  double xmax = -std::numeric_limits<double>::max();
  for (idx_t i : *this) {
    xmin = std::min(xmin, xy(i, XX) );
    xmax = std::max(xmax, xy(i, XX) );
  }
  comm.allReduceInPlace(xmin,eckit::mpi::min());
  comm.allReduceInPlace(xmax,eckit::mpi::max());

  size_t count = mesh_.nodes().size();
  size_t count_all = count;
  comm.allReduceInPlace(count_all, eckit::mpi::sum());


  for (int r = 0; r < mpi_size; ++r) {
      if (mpi_rank == r) {
          std::ofstream f(filepath.asString().c_str(), mpi_rank == 0? std::ios::trunc : std::ios::app);

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
          }
          f << "\n" "verts_" << r << " = [";
          for (idx_t i : static_cast<const container_t&>(*this)) {
              f << "\n  (" << xy(i, XX) << ", " << xy(i, YY) << "), ";
          }
          f << "\n]"
               "\n" ""
               "\n" "codes_" << r << " = [Path.MOVETO]"
               "\n" "codes_" << r << ".extend([Path.LINETO] * " << (size()-2) << ")"
               "\n" "codes_" << r << ".extend([Path.CLOSEPOLY])"
               "\n" ""
               "\n" "count_" << r << " = " << count <<
               "\n" "count_all_" << r << " = " << count_all <<
               "\n" ""
               "\n" "x_" << r << " = ["; for (size_t i=0; i<count; ++i) { f << xy(i, XX) << ", "; } f << "]"
               "\n" "y_" << r << " = ["; for (size_t i=0; i<count; ++i) { f << xy(i, YY) << ", "; } f << "]"
               "\n"
               "\n" "c = cycol()"
               "\n" "ax.add_patch(patches.PathPatch(Path(verts_" << r << ", codes_" << r << "), facecolor=c, color=c, alpha=0.3, lw=1))"
               "\n" "ax.scatter(x_" << r << ", y_" << r << ", color=c, marker='o')"
               "\n" "";
          if (mpi_rank == mpi_size - 1) {
              f << "\n" "ax.set_xlim(  0-5, 360+5)"  // ( "<<xmin<<"-5, "<<xmax<<"+5)
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


bool PartitionPolygon::containsPointInLonLatGeometry(const std::vector<PointLonLat>& points, const PointLonLat& P) const {
    ASSERT(points.size());

    // winding number
    int wn = 0;

    // loop on polygon edges
    for (size_t i = 1; i < points.size(); ++i) {
        const PointLonLat& A = points[i-1];
        const PointLonLat& B = points[ i ];

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


bool PartitionPolygon::containsPointInSphericalGeometry(const std::vector<PointLonLat>& points, const PointLonLat& P) const {
    ASSERT(points.size());

    // winding number
    int wn = 0;

    // loop on polygon edges
    for (size_t i = 1; i < points.size(); ++i) {
        const PointLonLat& A = points[i-1];
        const PointLonLat& B = points[ i ];

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


}  // namespace mesh
}  // namespace atlas

