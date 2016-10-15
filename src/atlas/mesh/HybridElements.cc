/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include "atlas/internals/atlas_config.h"
#include "atlas/internals/atlas_defines.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/field/Field.h"
#include "atlas/array/IndexView.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/ErrorHandling.h"
#include "eckit/log/Bytes.h"
#include "atlas/runtime/Log.h"

#ifdef ATLAS_HAVE_FORTRAN
#define FORTRAN_BASE 1
#define TO_FORTRAN +1
#else
#define FORTRAN_BASE 0
#endif

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

HybridElements::HybridElements() :
  size_(0),
  elements_size_(),
  elements_begin_(1,0ul),
  type_idx_()
{
  global_index_ = &add( field::Field::create<gidx_t>("glb_idx",   array::make_shape(size())) );
  remote_index_ = &add( field::Field::create<int   >("remote_idx",array::make_shape(size())) );
  partition_    = &add( field::Field::create<int   >("partition", array::make_shape(size())) );
  halo_         = &add( field::Field::create<int   >("halo",      array::make_shape(size())) );

  node_connectivity_ = &add( new Connectivity("node") );
  edge_connectivity_ = &add( new Connectivity("edge") );
  cell_connectivity_ = &add( new Connectivity("cell") );
}

HybridElements::~HybridElements()
{
}

field::Field& HybridElements::add( field::Field* field )
{
  ASSERT( field != NULL );
  ASSERT( ! field->name().empty() );

  if( has_field(field->name()) ) {
    std::stringstream msg;
    msg << "Trying to add field '"<<field->name()<<"' to HybridElements, but HybridElements already has a field with this name.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_[field->name()] = eckit::SharedPtr<field::Field>(field);
  return *field;
}

void HybridElements::resize( size_t size )
{
  size_ = size;
  for( FieldMap::iterator it = fields_.begin(); it != fields_.end(); ++it )
  {
    field::Field& field = *it->second;
    array::ArrayShape shape = field.shape();
    shape[0] = size_;
    field.resize(shape);
  }
}

void HybridElements::remove_field(const std::string& name)
{
  if( ! has_field(name) )
  {
    std::stringstream msg;
    msg << "Trying to remove field `"<<name<<"' in Nodes, but no field with this name is present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_.erase(name);
}


const field::Field& HybridElements::field(const std::string& name) const
{
  if( ! has_field(name) )
  {
    std::stringstream msg;
    msg << "Trying to access field `"<<name<<"' in Nodes, but no field with this name is present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  return *fields_.find(name)->second;
}

field::Field& HybridElements::field(const std::string& name)
{
  return const_cast<field::Field&>(static_cast<const HybridElements*>(this)->field(name));
}

const field::Field& HybridElements::field(size_t idx) const
{
  ASSERT(idx < nb_fields());
  size_t c(0);
  for( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it )
  {
    if( idx == c )
    {
      const field::Field& field = *it->second;
      return field;
    }
    c++;
  }
  throw eckit::SeriousBug("Should not be here!",Here());
}

field::Field& HybridElements::field(size_t idx)
{
  return const_cast<field::Field&>(static_cast<const HybridElements*>(this)->field(idx));
}


HybridElements::Connectivity& HybridElements::add( Connectivity *connectivity )
{
  connectivities_[connectivity->name()] = connectivity;
  return *connectivity;
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

  type_idx_.resize(new_size);

  for( size_t e=old_size; e<new_size; ++e ) {
    type_idx_[e] = element_types_.size();
  }

  elements_begin_.push_back(new_size);
  elements_size_.push_back(nb_elements);

  element_types_.push_back( etype );
  elements_.resize(element_types_.size());
  for( size_t t=0; t<nb_types(); ++t ) {
    if( elements_[t] ) elements_[t]->rebuild();
    else               elements_[t].reset( new Elements(*this,t) );
  }

//  for( size_t t=0; t<nb_types()-1; ++t )
//  {
//    elements_[t]->rebuild();
//  }
//  element_types_.push_back( etype );
//  elements_.push_back( eckit::SharedPtr<Elements>(new Elements(*this,type_idx_.back())) );


  node_connectivity_->add(nb_elements,nb_nodes,connectivity,fortran_array);

  resize( new_size );
  return element_types_.size()-1;
}

size_t HybridElements::add( const ElementType* element_type, size_t nb_elements )
{
  eckit::SharedPtr<const ElementType> etype ( element_type );

  size_t old_size=size();
  size_t new_size = old_size+nb_elements;

  size_t nb_nodes = etype->nb_nodes();

  type_idx_.resize(new_size);

  for( size_t e=old_size; e<new_size; ++e ) {
    type_idx_[e] = element_types_.size();
  }

  elements_begin_.push_back(new_size);
  elements_size_.push_back(nb_elements);

  element_types_.push_back( etype );
  elements_.resize(element_types_.size());
  for( size_t t=0; t<nb_types(); ++t ) {
    elements_[t].reset( new Elements(*this,t) );
  }

  node_connectivity_->add(nb_elements,nb_nodes);
  resize( new_size );
  return element_types_.size()-1;
}

size_t HybridElements::add( const Elements& elems )
{
  bool fortran_array = true;
  return add( &elems.element_type(), elems.size(), elems.node_connectivity().data(), fortran_array );
}


const std::string& HybridElements::name(size_t elem_idx) const { return element_types_[type_idx_[elem_idx]]->name(); }

size_t HybridElements::elemtype_nb_nodes(size_t elem_idx) const
{
  return element_type( type_idx(elem_idx) ).nb_nodes();
}

size_t HybridElements::elemtype_nb_edges(size_t elem_idx) const
{
  return element_type( type_idx(elem_idx) ).nb_edges();
}

void HybridElements::insert( size_t type_idx, size_t position, size_t nb_elements )
{
  type_idx_.insert(type_idx_.begin()+position,nb_elements,type_idx);
  elements_size_[type_idx] += nb_elements;
  for( size_t jtype=type_idx+1; jtype<nb_types()+1; ++jtype )
    elements_begin_[jtype] += nb_elements;
  for( size_t t=0; t<nb_types(); ++t ) {
    elements_[t]->rebuild();
  }
  node_connectivity_->insert(position,nb_elements,element_types_[type_idx]->nb_nodes());

  size_+=nb_elements;
  for( FieldMap::iterator it = fields_.begin(); it != fields_.end(); ++it )
  {
    field::Field& field = *it->second;
    field.insert(position,nb_elements);
  }

}


//-----------------------------------------------------------------------------

void HybridElements::clear()
{
  resize(0);
  for( ConnectivityMap::iterator it = connectivities_.begin(); it!=connectivities_.end(); ++it )
  {
    it->second->clear();
  }
  size_=0;
  elements_size_.clear();
  elements_begin_.resize(1);
  elements_begin_[0]=0;
  element_types_.clear();
  type_idx_.clear();
  elements_.clear();
}


size_t HybridElements::footprint() const {
  size_t size = sizeof(*this);
  for( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it ) {
    size += (*it).second->footprint();
  }
  for( ConnectivityMap::const_iterator it = connectivities_.begin(); it != connectivities_.end(); ++it ) {
    size += (*it).second->footprint();
  }
  size += elements_size_.capacity() * sizeof(size_t);
  size += elements_begin_.capacity() * sizeof(size_t);
  
  size += metadata_.footprint();
  
  return size;
}

//-----------------------------------------------------------------------------

extern "C" {
HybridElements* atlas__mesh__HybridElements__create()
{
  HybridElements* This = 0;
  ATLAS_ERROR_HANDLING( This = new HybridElements() );
  return This;
}

void atlas__mesh__HybridElements__delete(HybridElements* This)
{
  ATLAS_ERROR_HANDLING( delete This );
}

MultiBlockConnectivity* atlas__mesh__HybridElements__node_connectivity(HybridElements* This)
{
  MultiBlockConnectivity* connectivity(0);
  ATLAS_ERROR_HANDLING( connectivity = &This->node_connectivity() );
  return connectivity;
}

MultiBlockConnectivity* atlas__mesh__HybridElements__edge_connectivity(HybridElements* This)
{
  MultiBlockConnectivity* connectivity(0);
  ATLAS_ERROR_HANDLING( connectivity = &This->edge_connectivity() );
  return connectivity;
}

MultiBlockConnectivity* atlas__mesh__HybridElements__cell_connectivity(HybridElements* This)
{
  MultiBlockConnectivity* connectivity(0);
  ATLAS_ERROR_HANDLING( connectivity = &This->cell_connectivity() );
  return connectivity;
}

size_t atlas__mesh__HybridElements__size(const HybridElements* This)
{
  return This->size();
}

void atlas__mesh__HybridElements__add_elements(HybridElements* This, ElementType* elementtype, size_t nb_elements)
{
  This->add(elementtype,nb_elements);
}

void atlas__mesh__HybridElements__add_elements_with_nodes(HybridElements*This, ElementType* elementtype, size_t nb_elements, int node_connectivity[], int fortran_array)
{
  This->add(elementtype,nb_elements,node_connectivity,fortran_array);
}

int atlas__mesh__HybridElements__has_field(const HybridElements* This, char* name)
{
  ATLAS_ERROR_HANDLING( ASSERT(This!=0) );
  return This->has_field(std::string(name));
}

int atlas__mesh__HybridElements__nb_fields(const HybridElements* This)
{
  ATLAS_ERROR_HANDLING( ASSERT(This!=0) );
  return This->nb_fields();
}

int atlas__mesh__HybridElements__nb_types(const HybridElements* This)
{
  ATLAS_ERROR_HANDLING( ASSERT(This!=0) );
  return This->nb_types();
}

field::Field* atlas__mesh__HybridElements__field_by_idx(HybridElements* This, size_t idx)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->field(idx);
  );
  return field;
}

field::Field* atlas__mesh__HybridElements__field_by_name(HybridElements* This, char* name)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->field(std::string(name));
  );
  return field;
}

field::Field* atlas__mesh__HybridElements__global_index(HybridElements* This)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->global_index();
  );
  return field;

}

field::Field* atlas__mesh__HybridElements__remote_index(HybridElements* This)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->remote_index();
  );
  return field;
}

field::Field* atlas__mesh__HybridElements__partition(HybridElements* This)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->partition();
  );
  return field;
}

field::Field* atlas__mesh__HybridElements__halo(HybridElements* This)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->halo();
  );
  return field;
}

Elements* atlas__mesh__HybridElements__elements(HybridElements* This, size_t idx)
{
  Elements* elements(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    elements = &This->elements(idx);
  );
  return elements;
}

void atlas__mesh__HybridElements__add_field(HybridElements*This, field::Field* field)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    This->add(field);
  );
}


}

//-----------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas

#undef FORTRAN_BASE



