/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/array/MakeView.h"
#include "atlas/mesh/Polygon.h"
#include "atlas/field/Field.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/mesh.h"

using atlas::array::make_datatype;
using atlas::array::make_shape;

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

struct Edge : std::pair< size_t, size_t > {
    Edge(size_t _A, size_t _B) : std::pair< size_t, size_t >(_A, _B) {}
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


Polygon::Polygon(const detail::MeshImpl &mesh, size_t halo) :
  mesh_(mesh),
  halo_(halo) {
  const eckit::mpi::Comm& comm = parallel::mpi::comm();
  const int mpi_size = int(comm.size());

  // Polygon indices
  // - extract partition boundary edges by always attempting first to
  //   remove a reversed edge from a neighbouring element (if it exists)
  // - polygon can have multiple cycles, but must be a connected graph
  Poly poly;

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
                  Edge edge(size_t(conn(j, k)), size_t(conn(j, (k+1) % nb_nodes)));
                  if (!edges.erase(edge.reverse())) {
                      edges.insert(edge);
                  }
              }
          }
      }
  }

  while (!edges.empty()) {
      poly += Poly(edges);
  }
  ASSERT(poly.size());
  node_list_.assign(poly.begin(),poly.end());
}

size_t Polygon::footprint() const
{
    size_t size = sizeof(*this);
    size += node_list_.capacity()*sizeof(idx_t);
    return size;
}

void Polygon::print(std::ostream & out) const
{
    out << "polygon:{halo:"<<halo_<<",size:"<<size()<<",nodes:";
    for(size_t j=0; j<size()-1; ++j) {
      out << node_list_[j] << ",";
    }
    out << node_list_.back() << "}";
}

void Polygon::outputPythonScript(const eckit::PathName& filepath) const {
  const eckit::mpi::Comm& comm = atlas::parallel::mpi::comm();
  int mpi_rank = comm.rank();
  int mpi_size = comm.size();

  auto xy = array::make_view<double,2>( mesh_.nodes().xy() );

  double xmin =  std::numeric_limits<double>::max();
  double xmax = -std::numeric_limits<double>::max();
  for (size_t i = 0; i < size(); ++i) {
    xmin = std::min(xmin, xy(node_list_[i],XX) );
    xmax = std::max(xmax, xy(node_list_[i],XX) );
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
          for (size_t i = 0; i < size(); ++i) { f << "\n  (" << xy(node_list_[i],XX) << ", " << xy(node_list_[i],YY) << "), "; }
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

