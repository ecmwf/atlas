/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <stdexcept>
#include "eckit/exception/Exceptions.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/projection/Projection.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/field/Field.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/array/MakeView.h"

using atlas::grid::Grid;
using atlas::grid::Projection;

namespace atlas {
namespace mesh {

//----------------------------------------------------------------------------------------------------------------------

Mesh::Mesh( const Mesh& other ) :
    mesh_(other.mesh_) {
}
Mesh::Mesh( MeshImpl* mesh ) :
    mesh_(mesh) {
}
Mesh::Mesh( eckit::Stream& stream ) :
    mesh_( new MeshImpl(stream) ) {
}
Mesh::Mesh() :
    mesh_( new MeshImpl() ) {
}


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


mesh::Nodes& MeshImpl::createNodes(const Grid& grid)
{
  size_t nb_nodes = grid.size();
  nodes().resize(nb_nodes);

  array::ArrayView<double,2> lonlat = array::make_view<double,2>( nodes().lonlat() );
  array::ArrayView<double,2> geolonlat = array::make_view<double,2>( nodes().geolonlat() );
  size_t jnode(0);
  Projection projection = grid.projection();
  PointLonLat Pll;
  for( PointXY Pxy : grid ) {
    lonlat(jnode,0) = Pxy.x();
    lonlat(jnode,1) = Pxy.y();
    Pll = projection.lonlat(Pxy);
    geolonlat(jnode,0) = Pll.lon();
    geolonlat(jnode,1) = Pll.lat();
    ++jnode;
  }
  return nodes();
}

void MeshImpl::prettyPrint(std::ostream& os) const
{
}

void MeshImpl::print(std::ostream& os) const
{
}

size_t MeshImpl::footprint() const {
  size_t size = sizeof(*this);

  size += metadata_.footprint();
  if(nodes_)  size += nodes_->footprint();
  if(cells_)  size += cells_->footprint();
  if(facets_) size += facets_->footprint();
  if(ridges_) size += ridges_->footprint();
  if(peaks_)  size += peaks_->footprint();

  return size;
}


void MeshImpl::createElements()
{
  cells_.reset( new mesh::HybridElements() );
  facets_.reset( new mesh::HybridElements() );
  ridges_.reset( new mesh::HybridElements() );
  peaks_.reset( new mesh::HybridElements() );
  if( dimensionality_ == 2 )
    edges_ = facets_;
  else if( dimensionality_ == 3)
    edges_ = ridges_;
  else
    throw eckit::Exception("Invalid MeshImpl dimensionality",Here());

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

// C wrapper interfaces to C++ routines

MeshImpl* atlas__Mesh__new () {
	return new MeshImpl();
}

void atlas__Mesh__delete (MeshImpl* This) {
	delete This;
}

mesh::Nodes* atlas__Mesh__create_nodes (MeshImpl* This, int nb_nodes)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    This->nodes().resize(nb_nodes);
    return &This->nodes();
  );
  return NULL;
}

mesh::Nodes* atlas__Mesh__nodes (MeshImpl* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return &This->nodes();
  );
  return NULL;
}

mesh::Edges* atlas__Mesh__edges (MeshImpl* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return &This->edges();
  );
  return NULL;
}

mesh::Cells* atlas__Mesh__cells (MeshImpl* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return &This->cells();
  );
  return NULL;
}

size_t atlas__Mesh__footprint (MeshImpl* This) {
  size_t size(0);
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    size = This->footprint();
  );
  return size;
}



//----------------------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas
