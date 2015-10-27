/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/Parameters.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/Array.h"
#include "atlas/util/IndexView.h"

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

Connectivity::Connectivity(const Elements &elements, const int *values)
  : elements_(&elements), values_(values)
{
}

Elements::Elements(const Nodes& nodes) :
  nodes_(nodes),
  size_(0),
  nb_elements_(),
  type_begin_(),
  type_end_(),
  node_connectivity_( new ArrayT<int>() ),
  nb_nodes_( new ArrayT<int>() ),
  nb_edges_( new ArrayT<int>() ),
  type_idx_( new ArrayT<int>() ),
  node_connectivity_access_(*this,node_connectivity_->data())
{
}

Elements::~Elements() {}

//void Elements::add( ElementType* element_type, size_t nb_elements )
//{
//  eckit::SharedPtr<ElementType> etype ( element_type );

//  size_t nb_nodes = etype->nb_nodes();

//  size_t c(connectivity_->size());
//  connectivity_->resize(c+nb_elements*nb_nodes);

//  size_t n=nb_nodes_->size();
//  nb_nodes_->resize(n+nb_elements);

//  int* data_nb_nodes = nb_nodes_->data<int>();
//  for( size_t j=n; j<nb_nodes_->size(); ++j )
//    data_nb_nodes[j] = nb_nodes;

//  element_types_.push_back( etype );
//}

void Elements::add( ElementType* element_type, size_t nb_elements, const size_t connectivity[] )
{
  eckit::SharedPtr<ElementType> etype ( element_type );

  size_t nb_nodes = etype->nb_nodes();
  size_t nb_edges = etype->nb_edges();

  size_t c(node_connectivity_->size());
  node_connectivity_->resize(c+nb_elements*nb_nodes);

  size_t n=size();
  size_t new_size = n+nb_elements;
  nb_nodes_->resize(new_size);
  nb_edges_->resize(new_size);
  type_idx_->resize(new_size);

  int* data_nb_nodes = nb_nodes_->data();
  int* data_nb_edges = nb_edges_->data();
  int* data_type_idx = type_idx_->data();
  for( size_t j=n; j<new_size; ++j ) {
    data_nb_nodes[j] = nb_nodes;
    data_nb_edges[j] = nb_edges;
    data_type_idx[j] = element_types_.size();
  }

  int* conn = node_connectivity_->data();
  for(size_t j=0; c<node_connectivity_->size(); ++c, ++j) {
    conn[c] = connectivity[j];
  }

  type_begin_.push_back(size_);
  element_begin_.reserve(element_begin_.size()+nb_elements*nb_nodes);
  element_end_.  reserve(element_end_.  size()+nb_elements*nb_nodes);
  size_t element_size = size_ ? element_end_.back() : 0;
  for( size_t e=0; e<nb_elements; ++e ) {
    element_begin_.push_back(element_size+e*nb_nodes);
    element_end_.push_back(element_begin_.back()+nb_nodes);
  }
  size_ += nb_elements;
  type_end_.push_back(size_);
  nb_elements_.push_back(nb_elements);

  element_types_.push_back( etype );
  element_type_connectivity_.resize(element_types_.size());
  for( size_t t=0; t<nb_types(); ++t ) {
    element_type_connectivity_[t] = ElementTypeConnectivity(*this,t,node_connectivity_->data());
  }
  node_connectivity_access_ = Connectivity(*this,node_connectivity_->data());
}

const std::string& Elements::name(size_t elem_idx) const { return element_type( type_idx_->data()[elem_idx] ).name(); }
size_t Elements::nb_nodes(size_t elem_idx) const { return nb_nodes_->data()[elem_idx]; }
size_t Elements::nb_edges(size_t elem_idx) const { return nb_edges_->data()[elem_idx]; }

//-----------------------------------------------------------------------------

ElementTypeElements::ElementTypeElements(const Elements& elements, size_t type_idx)
  : elements_(&elements), type_idx_(type_idx)
{
  nb_nodes_ = elements.nb_nodes_->data()+elements.type_begin(type_idx_);
  nb_edges_ = elements.nb_edges_->data()+elements.type_begin(type_idx_);
}

const std::string& ElementTypeElements::name() const { return elements_->element_type(type_idx_).name(); }


extern "C" {

}

}  // namespace mesh
}  // namespace atlas

