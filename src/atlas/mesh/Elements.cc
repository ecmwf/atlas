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
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/Array.h"
#include "atlas/util/IndexView.h"
#include "atlas/mesh/Elements.h"


#ifdef ATLAS_HAVE_FORTRAN
#define FORTRAN_BASE 1
#define TO_FORTRAN +1
#else
#define FORTRAN_BASE 0
#define TO_FORTRAN
#endif

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

HybridElements::HybridElements() :
  size_(0),
  elements_size_(),
  elements_begin_(),
  node_connectivity_( new ArrayT<idx_t>() ),
  nodes_begin_(),
  nb_nodes_(),
  nb_edges_(),
  type_idx_(),
  node_connectivity_access_(node_connectivity_->data(),nodes_begin_.data())
{
}

HybridElements::~HybridElements()
{
  eckit::Log::info() << "destroying HybridElements" << std::endl;
}

size_t HybridElements::add( const ElementType* element_type, size_t nb_elements, const std::vector<idx_t> &connectivity )
{
  return add(element_type,nb_elements,connectivity.data());
}

size_t HybridElements::add( const ElementType* element_type, size_t nb_elements, const idx_t connectivity[] )
{
  return add(element_type,nb_elements,connectivity,false);
}

size_t HybridElements::add( const ElementType* element_type, size_t nb_elements, const idx_t connectivity[], bool fortran_array )
{
  eckit::SharedPtr<const ElementType> etype ( element_type );

  size_t old_size=size();
  size_t new_size = old_size+nb_elements;

  size_t nb_nodes = etype->nb_nodes();
  size_t nb_edges = etype->nb_edges();

  size_t c(node_connectivity_->size());
  node_connectivity_->resize(c+nb_elements*nb_nodes);

  nb_nodes_.resize(new_size);
  nb_edges_.resize(new_size);
  type_idx_.resize(new_size);
  nodes_begin_.resize(new_size);

  size_t nodes_end = old_size ? nodes_begin_[old_size-1]+nb_nodes : 0;
  for( size_t e=old_size; e<new_size; ++e ) {
    nb_nodes_[e] = nb_nodes;
    nb_edges_[e] = nb_edges;
    type_idx_[e] = element_types_.size();
    nodes_begin_[e] = nodes_end;
    nodes_end += nb_nodes;
  }

  idx_t *conn = node_connectivity_->data();
  idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
  for(size_t j=0; c<node_connectivity_->size(); ++c, ++j) {
    conn[c] = connectivity[j] + add_base ;
  }

  elements_begin_.push_back(old_size);
  elements_size_.push_back(nb_elements);

  element_types_.push_back( etype );
  element_type_connectivity_.resize(element_types_.size());
  elements_.resize(element_types_.size());
  for( size_t t=0; t<nb_types(); ++t ) {
    element_type_connectivity_[t] = Elements::Connectivity(
          node_connectivity_->data()+nodes_begin_[elements_begin_[t]],
          element_types_[t]->nb_nodes()
        );
    elements_[t] = Elements(*this,t);
  }
  node_connectivity_access_ = Connectivity(node_connectivity_->data(),nodes_begin_.data());
  return element_types_.size()-1;
}

size_t HybridElements::add( const Elements& elems )
{
  return add( &elems.element_type(), elems.size(), elems.hybrid_elements().node_connectivity_->data(), true );
}


const std::string& HybridElements::name(size_t elem_idx) const { return element_types_[type_idx_[elem_idx]]->name(); }
size_t HybridElements::nb_nodes(size_t elem_idx) const { return nb_nodes_[elem_idx]; }
size_t HybridElements::nb_edges(size_t elem_idx) const { return nb_edges_[elem_idx]; }

void HybridElements::set_node_connectivity( size_t elem_idx, const idx_t node_connectivity[] )
{
  idx_t *element_nodes = node_connectivity_->data()+nodes_begin_[elem_idx];
  for( size_t n=0; n<nb_nodes_[elem_idx]; ++n ) {
    element_nodes[n] = node_connectivity[n] TO_FORTRAN;
  }
}


void HybridElements::set_node_connectivity( size_t type_idx, size_t elem_idx, const idx_t node_connectivity[] )
{
  element_type_connectivity_[type_idx].set(elem_idx,node_connectivity);
}


//-----------------------------------------------------------------------------

Elements::Elements() : hybrid_elements_(0), type_idx_(0), nb_nodes_(0), nb_edges_(0), owns_elements_(false) {}

Elements::Elements(HybridElements &elements, size_t type_idx)
  : hybrid_elements_(&elements), type_idx_(type_idx), owns_elements_(false)
{
  nb_nodes_ = hybrid_elements_->element_type(type_idx_).nb_nodes();
  nb_edges_ = hybrid_elements_->element_type(type_idx_).nb_edges();
}

Elements::Elements(ElementType* element_type, size_t nb_elements, const std::vector<idx_t> &node_connectivity )
  : owns_elements_(true)
{
  hybrid_elements_ = new HybridElements();
  type_idx_ = hybrid_elements_->add(element_type,nb_elements,node_connectivity.data());
  nb_nodes_ = hybrid_elements_->element_type(type_idx_).nb_nodes();
  nb_edges_ = hybrid_elements_->element_type(type_idx_).nb_edges();
}

Elements::Elements( ElementType* element_type, size_t nb_elements, const idx_t node_connectivity[], bool fortran_array )
  : owns_elements_(true)
{
  hybrid_elements_ = new HybridElements();
  type_idx_ = hybrid_elements_->add(element_type,nb_elements,node_connectivity,fortran_array);
  nb_nodes_ = hybrid_elements_->element_type(type_idx_).nb_nodes();
  nb_edges_ = hybrid_elements_->element_type(type_idx_).nb_edges();
}


Elements::~Elements()
{
  if( owns_elements_ ) delete hybrid_elements_;
}

const std::string& Elements::name() const { return hybrid_elements_->element_type(type_idx_).name(); }

void Elements::set_node_connectivity( size_t elem_idx, const idx_t node_connectivity[] )
{
  return hybrid_elements_->set_node_connectivity(type_idx_,elem_idx,node_connectivity);
}

//-----------------------------------------------------------------------------

HybridConnectivity::HybridConnectivity(idx_t *values, size_t *offset)
  : values_(values), offset_(offset)
{
}

extern "C" {

}

}  // namespace mesh
}  // namespace atlas

