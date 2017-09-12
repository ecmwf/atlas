/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/Polygon.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

Polygon::Polygon(const detail::MeshImpl &mesh, size_t halo) :
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

  ASSERT(operator+=(Poly(edges)));
}

size_t Polygon::footprint() const
{
    size_t size = sizeof(*this);
    size += capacity()*sizeof(idx_t);
    return size;
}

void Polygon::print(std::ostream & out) const
{
    out << "polygon:{"
        <<  "halo:" << halo_
        << ",size:" << size()
        << ",nodes:" << static_cast<Poly>(*this)
        << "}";
}

void Polygon::outputPythonScript(const eckit::PathName& filepath) const {
  const eckit::mpi::Comm& comm = atlas::parallel::mpi::comm();
  int mpi_rank = comm.rank();
  int mpi_size = comm.size();

  auto xy = array::make_view<double,2>( mesh_.nodes().xy() );

  double xmin =  std::numeric_limits<double>::max();
  double xmax = -std::numeric_limits<double>::max();
  for (idx_t i : *this) {
    xmin = std::min(xmin, xy(i,XX) );
    xmax = std::max(xmax, xy(i,XX) );
  }
  comm.allReduceInPlace(xmin,eckit::mpi::min());
  comm.allReduceInPlace(xmax,eckit::mpi::max());


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
          for (idx_t i : *this) {
              f << "\n  (" << xy(i,XX) << ", " << xy(i,YY) << "), ";
          }
          f << "\n]"
               "\n" ""
               "\n" "codes_" << r << " = [Path.MOVETO]"
               "\n" "codes_" << r << ".extend([Path.LINETO] * " << (size()-2) << ")"
               "\n" "codes_" << r << ".extend([Path.CLOSEPOLY])"
               "\n" ""
               "\n" "c = cycol()"
               "\n" "ax.add_patch(patches.PathPatch(Path(verts_" << r << ", codes_" << r << "), facecolor=c, color=c, alpha=0.3, lw=1))"
               "\n" "";
          if (mpi_rank == mpi_size - 1) {
              f <<
                   "\n" "ax.set_xlim( "<<xmin<<"-5, "<<xmax<<"+5)"
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

}  // namespace mesh
}  // namespace atlas

