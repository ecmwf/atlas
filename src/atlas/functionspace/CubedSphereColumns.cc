/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field/Field.h"
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"

namespace atlas {
namespace functionspace {

// Helper functions to get fields.
namespace
{
template<typename BaseFunctionSpace>
Field getTij( const Mesh& mesh );

template<typename BaseFunctionSpace>
Field getGhost( const Mesh& mesh );

template<>
Field getTij<NodeColumns>( const Mesh& mesh ) {
  return mesh.nodes().field( "tij" );
}

template<>
Field getTij<CellColumns>( const Mesh& mesh ) {
  return mesh.cells().field( "tij" );
}

template<>
Field getGhost<NodeColumns>( const Mesh& mesh ) {
  return mesh.nodes().ghost();
}

template<>
Field getGhost<CellColumns>( const Mesh& mesh ) {
  return mesh.cells().halo();
}
}

// All constructors pass arguments through to BaseFunctionSpace, then construct
// CubedSphereColumnsImpl.
template<typename BaseFunctionSpace>
CubedSphereColumns<BaseFunctionSpace>::CubedSphereColumns() :
  BaseFunctionSpace (),
  cubedSphereColumnsHandle_( new detail::CubedSphereColumnsImpl() ) {}

template<typename BaseFunctionSpace>
CubedSphereColumns<BaseFunctionSpace>::CubedSphereColumns(
  const FunctionSpace& functionspace) :
  BaseFunctionSpace( functionspace ),
  cubedSphereColumnsHandle_( new detail::CubedSphereColumnsImpl(
    getTij<BaseFunctionSpace>( this->mesh() ),
    getGhost<BaseFunctionSpace>( this->mesh() ) ) ) {}

template<typename BaseFunctionSpace>
CubedSphereColumns<BaseFunctionSpace>::CubedSphereColumns(
  const Mesh& mesh, const eckit::Configuration& configuration) :
  BaseFunctionSpace( mesh, configuration ),
  cubedSphereColumnsHandle_( new detail::CubedSphereColumnsImpl(
    getTij<BaseFunctionSpace>( this->mesh() ),
    getGhost<BaseFunctionSpace>( this->mesh() ) ) ) {}

template<typename BaseFunctionSpace>
CubedSphereColumns<BaseFunctionSpace>::CubedSphereColumns(
  const Mesh& mesh) :
  BaseFunctionSpace( mesh ),
  cubedSphereColumnsHandle_( new detail::CubedSphereColumnsImpl(
    getTij<BaseFunctionSpace>( this->mesh() ),
    getGhost<BaseFunctionSpace>( this->mesh() ) ) ) {}

template<typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::i_begin(idx_t t) const {
  return cubedSphereColumnsHandle_.get()->i_begin(t);
}

template<typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::i_end(idx_t t) const {
  return cubedSphereColumnsHandle_.get()->i_end(t);
}

template<typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::j_begin(idx_t t) const {
  return cubedSphereColumnsHandle_.get()->j_begin(t);
}

template<typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::j_end(idx_t t) const {
  return cubedSphereColumnsHandle_.get()->j_end(t);
}

template<typename BaseFunctionSpace>
idx_t CubedSphereColumns<BaseFunctionSpace>::index(idx_t t, idx_t i, idx_t j) const {
  return cubedSphereColumnsHandle_.get()->index(t, i, j);
}

template<typename BaseFunctionSpace>
Field CubedSphereColumns<BaseFunctionSpace>::tij() const {
  return cubedSphereColumnsHandle_.get()->tij();
}

// Explicit instantiation of template classes.
template class CubedSphereColumns<CellColumns>;
template class CubedSphereColumns<NodeColumns>;

} // namespace functionspace
} // namespace atlas
