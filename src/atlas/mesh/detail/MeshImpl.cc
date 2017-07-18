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
  if(nodes_)           size += nodes_           ->footprint();
  if(cells_)           size += cells_           ->footprint();
  if(facets_)          size += facets_          ->footprint();
  if(ridges_)          size += ridges_          ->footprint();
  if(peaks_)           size += peaks_           ->footprint();
  if(partition_graph_) size += partition_graph_ ->footprint();
  for( const auto& polygon : polygons_ ) {
    if( polygon ) size += polygon->footprint();
  }

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

const PartitionGraph& MeshImpl::partitionGraph() const {
  if( not partition_graph_ ) {
    partition_graph_.reset( build_partition_graph(*this) );
  }
  return *partition_graph_;
}

PartitionGraph::Neighbours MeshImpl::nearestNeighbourPartitions() const {
  return partitionGraph().nearestNeighbours(partition());
}

const Polygon& MeshImpl::polygon(size_t halo) const {
  if( halo >= polygons_.size() ) {
    polygons_.resize(halo+1);
  }
  if( not polygons_[halo] ) {

    int mesh_halo = 0;
    metadata().get("halo",mesh_halo);
    if( halo > mesh_halo ) {
      throw eckit::Exception("Mesh does not contain a halo of size "+std::to_string(halo)+".", Here());
    }

    polygons_[halo].reset( new Polygon(*this, halo) );
  }
  return *polygons_[halo];
}


//----------------------------------------------------------------------------------------------------------------------


} // namespace detail
} // namespace mesh
} // namespace atlas
