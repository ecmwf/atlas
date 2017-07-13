/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iomanip>

#include "eckit/exception/Exceptions.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/log/Bytes.h"

#include "atlas/mesh/detail/MeshImpl.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Unique.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace mesh {
namespace detail {

//----------------------------------------------------------------------------------------------------------------------

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


PartitionGraph* build_partition_graph(const MeshImpl& mesh ) {

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
      auto patch = elements.view< int, 1 >(elements.field("patch"));

      const size_t nb_nodes = elements.nb_nodes();

      for (size_t j = 0; j < elements.size(); ++j) {
          for (size_t k = 0; k < nb_nodes && (patch(j) == 0); ++k) {
              Edge edge(size_t(conn(j, k)), size_t(conn(j, (k+1) % nb_nodes)));
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


  std::vector< double > polygon;
  polygon.reserve(poly.size()*2);

  auto xy = array::make_view< double, 2 >( mesh.nodes().xy() );

  for (size_t i = 0; i < poly.size(); ++i) {
      polygon.push_back( xy(poly[i],XX) );
      polygon.push_back( xy(poly[i],YY) );
  }
  ASSERT(polygon.size() >= 4);

  eckit::mpi::Buffer<double> recv_polygons(mpi_size);
  comm.allGatherv(polygon.begin(),polygon.end(),recv_polygons);

  using Polygon = std::vector<PointXY>;
  std::vector< Polygon > polygons(mpi_size);
  for( size_t p=0; p<mpi_size; ++p ) {
    for( size_t j=0; j<recv_polygons.counts[p]/2; ++j) {
      PointXY pxy(
        *(recv_polygons.begin()+recv_polygons.displs[p] + 2*j + XX),
        *(recv_polygons.begin()+recv_polygons.displs[p] + 2*j + YY)
      );
      polygons[p].push_back(pxy);
    }
  }

  std::map< uidx_t, std::set<size_t> > uid_2_parts;
  size_t jpart=0;
  for( const Polygon& _polygon : polygons ) {
    for( const PointXY& pxy : _polygon ) {
      PointLonLat pll = pxy;
      if( eckit::types::is_strictly_greater( 0., pll.lon() ) ) {
        pll.lon() += 360.;
      }
      if( eckit::types::is_approximately_greater_or_equal( pll.lon(), 360. ) ) {
        pll.lon() -= 360.;
      }
      uidx_t uid = util::unique_lonlat( pll.data() );
      uid_2_parts[uid].insert( jpart );
    }
    ++jpart;
  }
  std::vector< std::set<size_t> > graph(mpi_size);
  for( const auto& u2p : uid_2_parts ){
    const std::set<size_t>& parts = u2p.second;
    for( size_t jpart : parts ) {
      for( size_t ipart : parts ) {
        if( jpart != ipart ) {
          graph[jpart].insert(ipart);
        }
      }
    }
  }

  std::vector<size_t> counts(mpi_size);
  std::vector<size_t> displs(mpi_size);
  size_t values_size=0;
  for( size_t jpart=0; jpart<mpi_size; ++jpart ) {
    counts[jpart] = graph[jpart].size();
    displs[jpart] = values_size;
    values_size += counts[jpart];
  }
  std::vector<size_t> values; values.reserve(values_size);
  for( const std::set<size_t>& graph_node : graph ) {
    for( size_t v : graph_node ) {
      values.push_back(v);
    }
  }

  return new PartitionGraph( values.data(), mpi_size, displs.data(), counts.data() );
}


size_t PartitionGraph::footprint() const {
  size_t size = sizeof(*this);
  size += sizeof(size_t)*displs_.capacity();
  size += sizeof(size_t)*counts_.capacity();
  size += sizeof(size_t)*values_.capacity();
  return size;
}

size_t PartitionGraph::size() const {
  return displs_.size();
}

PartitionGraph::Neighbours PartitionGraph::nearestNeighbours(const size_t partition) const {
  return Neighbours(values_.data()+displs_[partition],values_.data()+displs_[partition]+counts_[partition]);
}

PartitionGraph::PartitionGraph() {
}

PartitionGraph::PartitionGraph( size_t values[], size_t rows, size_t displs[], size_t counts[] ) {
  displs_.assign(displs,displs+rows);
  counts_.assign(counts,counts+rows);
  values_.assign(values,values+displs[rows-1]+counts[rows-1]);
  maximum_nearest_neighbours_ = 0;
  for( size_t n : counts_ ) {
    maximum_nearest_neighbours_ = std::max(n,maximum_nearest_neighbours_);
  }
}

size_t PartitionGraph::maximumNearestNeighbours() const {
  return maximum_nearest_neighbours_;
}

void PartitionGraph::print(std::ostream& os) const
{
  for( size_t jpart=0; jpart<size(); ++jpart ) {
    Log::info() << std::setw(3) << jpart << " : ";
    for( size_t v : nearestNeighbours(jpart) ) {
      Log::info() << std::setw(3) << v << " ";
    }
    Log::info() << '\n';
  }
  Log::info() << "partition graph maximum neighbours = " << maximumNearestNeighbours() << '\n';
  Log::info() << "partition graph footprint = " << eckit::Bytes(footprint());
}

PartitionGraph::operator bool() const {
  return size();
}

std::ostream& operator<<(std::ostream& s, const PartitionGraph& p) {
    p.print(s);
    return s;
}

//----------------------------------------------------------------------------------------------------------------------


} // namespace detail
} // namespace mesh
} // namespace atlas
