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
  elements_begin_(1,0ul),
  node_connectivity_array_( new ArrayT<idx_t>() ),
  nodes_begin_(1,0ul),
  nb_nodes_(),
  nb_edges_(),
  type_idx_()
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

  size_t c(node_connectivity_array_->size());
  node_connectivity_array_->resize(c+nb_elements*nb_nodes);

  nb_nodes_.resize(new_size);
  nb_edges_.resize(new_size);
  type_idx_.resize(new_size);
  nodes_begin_.resize(new_size+1);

  for( size_t e=old_size; e<new_size; ++e ) {
    nb_nodes_[e] = nb_nodes;
    nb_edges_[e] = nb_edges;
    type_idx_[e] = element_types_.size();
    nodes_begin_[e+1] = nodes_begin_[e]+nb_nodes;
  }

  idx_t *conn = node_connectivity_array_->data();
  idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
  for(size_t j=0; c<node_connectivity_array_->size(); ++c, ++j) {
    conn[c] = connectivity[j] + add_base ;
  }

  elements_begin_.push_back(new_size);
  elements_size_.push_back(nb_elements);

  element_types_.push_back( etype );
  elements_.resize(element_types_.size());
  for( size_t t=0; t<nb_types(); ++t ) {
    elements_[t] = Elements(*this,t);
  }
  node_connectivity_.reset( new HybridElements::Connectivity(
                                     node_connectivity_array_->data(),
                                     size_,
                                     nodes_begin_.data(),
                                     element_types_.size(),
                                     elements_begin_.data())
                                   );
  return element_types_.size()-1;
}

size_t HybridElements::add( const Elements& elems )
{
  return add( &elems.element_type(), elems.size(), elems.hybrid_elements().node_connectivity_array_->data(), true );
}


const std::string& HybridElements::name(size_t elem_idx) const { return element_types_[type_idx_[elem_idx]]->name(); }
size_t HybridElements::nb_nodes(size_t elem_idx) const { return nb_nodes_[elem_idx]; }
size_t HybridElements::nb_edges(size_t elem_idx) const { return nb_edges_[elem_idx]; }

void HybridElements::set_node_connectivity( size_t elem_idx, const idx_t node_connectivity[] )
{
  idx_t *element_nodes = node_connectivity_array_->data()+nodes_begin_[elem_idx];
  for( size_t n=0; n<nb_nodes_[elem_idx]; ++n ) {
    element_nodes[n] = node_connectivity[n] TO_FORTRAN;
  }
}


void HybridElements::set_node_connectivity( size_t type_idx, size_t elem_idx, const idx_t node_connectivity[] )
{
  node_connectivity_.get()->block_connectivity(type_idx).set(elem_idx,node_connectivity);
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

HybridConnectivity::  HybridConnectivity( idx_t *values, size_t rows, size_t *row_offset, size_t blocks, size_t *block_offset )
  : owns_(false),
    values_(values),
    rows_(rows),
    row_offset_(row_offset),
    blocks_(blocks),
    block_offset_(block_offset)
{
  regenerate_block_connectivity();
}


HybridConnectivity::HybridConnectivity() :
  owns_(true),
  values_(0),
  rows_(0),
  row_offset_(0),
  blocks_(0),
  block_offset_(0),
  owned_row_offset_(1,0ul),
  owned_block_offset_(1,0ul)
{}

HybridConnectivity::~HybridConnectivity() {}

void HybridConnectivity::add(size_t rows, size_t cols, const idx_t values[], bool fortran_array )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  size_t old_size = owned_values_.size();
  size_t new_size = old_size + rows*cols;
  size_t new_rows = rows_+rows;
  owned_row_offset_.resize(new_rows+1);
  for(size_t j=0; rows_<new_rows; ++rows_, ++j) {
    owned_row_offset_[rows_+1] = owned_row_offset_[rows_]+cols;
  }

  owned_values_.resize(new_size);
  idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
  for(size_t j=0, c=old_size; c<new_size; ++c, ++j) {
    owned_values_[c] = values[j] + add_base;
  }

  owned_block_offset_.push_back(new_rows);
  blocks_++;

  values_ = owned_values_.data();
  row_offset_ = owned_row_offset_.data();
  block_offset_ = owned_block_offset_.data();

  regenerate_block_connectivity();
}

void HybridConnectivity::add( const BlockConnectivity& block )
{
  if( !owns_ ) throw eckit::AssertionFailed("HybridConnectivity must be owned to be resized directly");
  bool fortran_array = FORTRAN_BASE;
  add(block.rows(),block.cols(),block.values_,fortran_array);
}


void HybridConnectivity::regenerate_block_connectivity()
{
  block_connectivity_.resize(blocks_);
  for( size_t b=0; b<blocks_; ++b )
  {
    block_connectivity_[b].reset(
       new BlockConnectivity(
        block_offset_[b+1]-block_offset_[b], // rows
        row_offset_[block_offset_[b]+1]-row_offset_[block_offset_[b]],  // cols
        values_+row_offset_[block_offset_[b]]) );
  }
}

BlockConnectivity::BlockConnectivity() :
  owns_(true), values_(0), rows_(0), cols_(0)
{
}

void BlockConnectivity::add(size_t rows, size_t cols, const idx_t *values, bool fortran_array)
{
  if( !owns_ )
    throw eckit::AssertionFailed("BlockConnectivity must be owned to be resized directly");
  if( cols_!=0 && cols_!=cols)
    throw eckit::AssertionFailed("Cannot add values with different cols than already existing in BlockConnectivity");

  size_t old_size = rows_*cols_;
  size_t new_size = old_size+rows*cols;
  owned_values_.resize(new_size);
  idx_t add_base = fortran_array ? 0 : FORTRAN_BASE;
  for( size_t j=0, c=old_size; c<new_size; ++c, ++j ) {
    owned_values_[c] = values[j] + add_base;
  }

  values_=owned_values_.data();
  rows_+=rows;
  cols_=cols;
}


extern "C" {

}

}  // namespace mesh
}  // namespace atlas

