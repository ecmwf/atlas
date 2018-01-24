/*
 * (C) Copyright 2013 ECMWF.
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

PartitionGraph* build_partition_graph(const MeshImpl& mesh ) {

  const eckit::mpi::Comm& comm = parallel::mpi::comm();
  const int mpi_size = int(comm.size());

  const util::Polygon& poly = mesh.polygon();

  std::vector< double > polygon;
  polygon.reserve(poly.size()*2);

  auto xy = array::make_view< double, 2 >( mesh.nodes().xy() );

  for( idx_t node : poly ) {
      polygon.push_back( xy(node,XX) );
      polygon.push_back( xy(node,YY) );
  }
  ASSERT(polygon.size() >= 4);

  eckit::mpi::Buffer<double> recv_polygons(mpi_size);
  comm.allGatherv(polygon.begin(),polygon.end(),recv_polygons);

  using PolygonXY = std::vector<PointXY>;
  std::vector< PolygonXY > polygons(mpi_size);
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
  for( const PolygonXY& _polygon : polygons ) {
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

  for( size_t jpart=0; jpart<rows; ++jpart ) {
    for( size_t neighbour : nearestNeighbours(jpart) ) {
      bool found(false);
      for( size_t nextneighbour : nearestNeighbours(neighbour) ) {
        if( nextneighbour == jpart )
          found = true;
      }
      if( not found ) {
        values_.insert( values_.begin()+displs_[neighbour]+counts_[neighbour], jpart );
        counts_[neighbour]++;
        for( size_t j=neighbour+1; j<rows; ++j ) {
          displs_[j]++;
        }
      }
    }
  }

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
