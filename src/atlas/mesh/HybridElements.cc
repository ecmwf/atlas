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
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/Field.h"
#include "atlas/atlas_config.h"
#include "atlas/atlas_defines.h"

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
  global_index_ = &add( Field::create<gidx_t>("glb_idx",   make_shape(size(),1)) );
  remote_index_ = &add( Field::create<int   >("remote_idx",make_shape(size(),1)) );
  partition_    = &add( Field::create<int   >("partition", make_shape(size(),1)) );
  halo_         = &add( Field::create<int   >("halo",     make_shape(size(),1)) );
  
  node_connectivity_ = &add( "node", new Connectivity() );  
  edge_connectivity_ = &add( "edge", new Connectivity() );  
  cell_connectivity_ = &add( "cell", new Connectivity() );  
}

HybridElements::~HybridElements()
{
}

Field& HybridElements::add( Field* field )
{
  ASSERT( field != NULL );
  ASSERT( ! field->name().empty() );

  if( has_field(field->name()) ) {
    std::stringstream msg;
    msg << "Trying to add field '"<<field->name()<<"' to Nodes, but Nodes already has a field with this name.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_[field->name()] = eckit::SharedPtr<Field>(field);
  return *field;
}

void HybridElements::resize( size_t size )
{
  size_ = size;
  for( FieldMap::iterator it = fields_.begin(); it != fields_.end(); ++it )
  {
    Field& field = *it->second;
    ArrayShape shape = field.shape();
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


const Field& HybridElements::field(const std::string& name) const
{
  if( ! has_field(name) )
  {
    std::stringstream msg;
    msg << "Trying to access field `"<<name<<"' in Nodes, but no field with this name is present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  return *fields_.find(name)->second;
}

Field& HybridElements::field(const std::string& name)
{
  return const_cast<Field&>(static_cast<const HybridElements*>(this)->field(name));
}

const Field& HybridElements::field(size_t idx) const
{
  ASSERT(idx < nb_fields());
  size_t c(0);
  for( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it )
  {
    if( idx == c )
    {
      const Field& field = *it->second;
      return field;
    }
    c++;
  }
  eckit::SeriousBug("Should not be here!",Here());
  static Field* ret;
  return *ret;
}

Field& HybridElements::field(size_t idx)
{
  return const_cast<Field&>(static_cast<const HybridElements*>(this)->field(idx));
}


HybridElements::Connectivity& HybridElements::add( const std::string& name, Connectivity *connectivity )
{
  connectivities_[name] = connectivity;
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
  size_t nb_edges = etype->nb_edges();

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
  
  node_connectivity_->add(nb_elements,nb_nodes,connectivity,fortran_array);
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

//-----------------------------------------------------------------------------

extern "C" {

}

//-----------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas

#undef FORTRAN_BASE



