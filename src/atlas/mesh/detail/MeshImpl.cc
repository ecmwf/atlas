/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/detail/MeshImpl.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Unique.h"
#include "atlas/runtime/Log.h"
#include "eckit/log/Bytes.h"
#include <iomanip>

using atlas::Grid;
using atlas::Projection;

namespace atlas {
namespace mesh {
namespace detail {

//----------------------------------------------------------------------------------------------------------------------

MeshImpl::MeshImpl(eckit::Stream& s)
{
    NOTIMP;
}

void MeshImpl::encode(eckit::Stream& s) const {
    NOTIMP;
}

MeshImpl::MeshImpl():
  dimensionality_(2)
{
  nodes_.reset( new mesh::Nodes() );
  createElements();
}


MeshImpl::~MeshImpl()
{
}

void MeshImpl::print(std::ostream& os) const
{
}

size_t MeshImpl::footprint() const {
  size_t size = sizeof(*this);

  size += metadata_.footprint();
  if(nodes_)  size += nodes_ ->footprint();
  if(cells_)  size += cells_ ->footprint();
  if(facets_) size += facets_->footprint();
  if(ridges_) size += ridges_->footprint();
  if(peaks_)  size += peaks_ ->footprint();

  return size;
}


void MeshImpl::createElements()
{
  cells_ .reset( new HybridElements() );
  facets_.reset( new HybridElements() );
  ridges_.reset( new HybridElements() );
  peaks_ .reset( new HybridElements() );
  if( dimensionality_ == 2 )
    edges_ = facets_;
  else if( dimensionality_ == 3)
    edges_ = ridges_;
  else
    throw eckit::Exception("Invalid Mesh dimensionality",Here());

  ASSERT( edges_.owners() == 2 );
}

bool MeshImpl::generated() const {
  return ! (cells_->size() == 0 && facets_->size() == 0 && ridges_->size() == 0 && peaks_->size() == 0);
}

void MeshImpl::setProjection(const Projection& projection) {
  projection_ = projection;
}

size_t MeshImpl::nb_partitions() const {
  return parallel::mpi::comm().size();
}

size_t MeshImpl::partition() const {
  return parallel::mpi::comm().rank();
}

void MeshImpl::cloneToDevice() const {
  if( nodes_  ) nodes_ ->cloneToDevice();
  if( cells_  ) cells_ ->cloneToDevice();
  if( facets_ ) facets_->cloneToDevice();
  if( ridges_ ) ridges_->cloneToDevice();
  if( peaks_  ) peaks_ ->cloneToDevice();
}

void MeshImpl::cloneFromDevice() const {
  if( nodes_  ) nodes_ ->cloneFromDevice();
  if( cells_  ) cells_ ->cloneFromDevice();
  if( facets_ ) facets_->cloneFromDevice();
  if( ridges_ ) ridges_->cloneFromDevice();
  if( peaks_  ) peaks_ ->cloneFromDevice();
}

void MeshImpl::syncHostDevice() const {
  if( nodes_  ) nodes_ ->syncHostDevice();
  if( cells_  ) cells_ ->syncHostDevice();
  if( facets_ ) facets_->syncHostDevice();
  if( ridges_ ) ridges_->syncHostDevice();
  if( peaks_  ) peaks_ ->syncHostDevice();
}

//----------------------------------------------------------------------------------------------------------------------

/// TODO: Move to other file

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

const PartitionGraph& MeshImpl::partitionGraph() const {
  if( not partition_graph_ ) {
    build_partition_graph();
  }
  return *partition_graph_;
}


void MeshImpl::build_partition_graph() const {
  
  const eckit::mpi::Comm& comm = parallel::mpi::comm();
  const int mpi_size = int(comm.size());

  // Polygon indices
  // - extract partition boundary edges by always attempting first to
  //   remove a reversed edge from a neighbouring element (if it exists)
  // - polygon can have multiple cycles, but must be a connected graph
  Poly poly;

  edge_set_t edges;
  for (size_t t = 0; t < cells().nb_types(); ++t) {
      const Elements& elements = cells().elements(t);

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

  auto xy = array::make_view< double, 2 >( nodes().xy() );

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

  partition_graph_.reset( new PartitionGraph( values.data(), mpi_size, displs.data(), counts.data() ) );
}


size_t PartitionGraph::footprint() const {
  return sizeof(displs_)+sizeof(size_t)*displs_.capacity()
       + sizeof(counts_)+sizeof(size_t)*counts_.capacity()
       + sizeof(values_)+sizeof(size_t)*values_.capacity()
       + sizeof(maximum_nearest_neighbours_);
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
  Log::info() << "maximum neighbours = " << maximumNearestNeighbours() << '\n';
  Log::info() << "footprint = " << eckit::Bytes(footprint());
}
PartitionGraph::operator bool() const {
  return size();
}

PartitionGraph::Neighbours MeshImpl::nearestNeighbourPartitions() const {
  return partitionGraph().nearestNeighbours(partition());
}



//----------------------------------------------------------------------------------------------------------------------


} // namespace detail
} // namespace mesh
} // namespace atlas
