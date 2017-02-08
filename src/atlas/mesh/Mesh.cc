/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/field/Field.h"
#include "atlas/internals/Parameters.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/ErrorHandling.h"

namespace atlas {
namespace mesh {

//----------------------------------------------------------------------------------------------------------------------

Mesh* Mesh::create( const eckit::Parametrisation& params )
{
  return new Mesh(params);
}



Mesh* Mesh::create( const grid::Grid& grid, const eckit::Parametrisation& params )
{
    return new Mesh(grid,params);
}

Mesh::Mesh(eckit::Stream& s)
{
    NOTIMP;
}

void Mesh::encode(eckit::Stream& s) const {
    NOTIMP;
}

Mesh::Mesh( const eckit::Parametrisation& ):
  dimensionality_(2)
{
  nodes_.reset( new mesh::Nodes() );
  createElements();
}

Mesh::Mesh(const grid::Grid& grid, const eckit::Parametrisation& ) :
   dimensionality_(2)
{
  nodes_.reset( new mesh::Nodes() );
  createNodes(grid);
  createElements();
}


Mesh::~Mesh()
{
}


mesh::Nodes& Mesh::createNodes(const grid::Grid& g)
{
  size_t nb_nodes = g.npts();
  nodes().resize(nb_nodes);
  g.fillLonLat(nodes().lonlat().data<double>(), nb_nodes*2);
  return nodes();
}

void Mesh::prettyPrint(std::ostream& os) const
{
}

void Mesh::print(std::ostream& os) const
{
}


void Mesh::createElements()
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
    throw eckit::Exception("Invalid Mesh dimensionality",Here());

  ASSERT( edges_.owners() == 2 );
}

bool Mesh::generated() const {
  return ! (cells_->size() == 0 && facets_->size() == 0 && ridges_->size() == 0 && peaks_->size() == 0);
}

void Mesh::setProjection(const grid::projection::Projection& prj) {
  projection_=prj.clone();
}

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines

Mesh* atlas__Mesh__new () {
	return new Mesh();
}

void atlas__Mesh__delete (Mesh* This) {
	delete This;
}

mesh::Nodes* atlas__Mesh__create_nodes (Mesh* This, int nb_nodes)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    This->nodes().resize(nb_nodes);
    return &This->nodes();
  );
  return NULL;
}

mesh::Nodes* atlas__Mesh__nodes (Mesh* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return &This->nodes();
  );
  return NULL;
}

mesh::Edges* atlas__Mesh__edges (Mesh* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return &This->edges();
  );
  return NULL;
}

mesh::Cells* atlas__Mesh__cells (Mesh* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return &This->cells();
  );
  return NULL;
}


//----------------------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas
