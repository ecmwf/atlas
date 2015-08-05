/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/Nodes.h"
#include "atlas/field/FieldT.h"

using atlas::field::FieldT;

namespace atlas {

//------------------------------------------------------------------------------------------------------

Nodes::Nodes(size_t _size) :
  FunctionSpace("nodes",make_shape(_size,FunctionSpace::UNDEF_VARS))
{
  global_index_ = &add( new FieldT<gidx_t>("glb_idx",   make_shape(size(),1)) );
  remote_index_ = &add( new FieldT<int   >("remote_idx",make_shape(size(),1)) );
  partition_    = &add( new FieldT<int   >("partition", make_shape(size(),1)) );
  ghost_        = &add( new FieldT<int   >("ghost",     make_shape(size(),1)) );
  halo_         = &add( new FieldT<int   >("halo",      make_shape(size(),1)) );
  topology_     = &add( new FieldT<int   >("topology",  make_shape(size(),1)) );
  lonlat_       = &add( new FieldT<double>("lonlat",    make_shape(size(),2)) );

  add( new FieldT<int>("flags", make_shape(size(),1)) );


  ArrayView<gidx_t,1> glb_idx( global_index() );
  ArrayView<int   ,1> part( partition() );
  ArrayView<int   ,1> flags( field("flags") );

  for(size_t n=0; n<size(); ++n)
  {
    glb_idx(n) = 1+n;
    part(n) = eckit::mpi::rank();
    flags(n) = 0;
  }
  metadata().set("nb_owned",size());
  metadata().set<long>("type",static_cast<int>(Entity::NODES));
}


Field& Nodes::add( Field* field )
{
  ASSERT( field != NULL );
  ASSERT( ! field->name().empty() );

  if( has_field(field->name()) ) {
    std::stringstream msg;
    msg << "Trying to add field '"<<field->name()<<"' to Nodes, but Nodes already has a field with this name.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_[field->name()] = eckit::SharedPtr<Field>(field);


  field->set_function_space(*this);


  return *field;
}

const Field& Nodes::field(const std::string& name) const
{
  if( ! has_field(name) )
  {
    std::stringstream msg;
    msg << "Trying to access field `"<<name<<"' in Nodes, but no field with this name is present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  return *fields_.find(name)->second;
}

Field& Nodes::field(const std::string& name)
{
  return const_cast<Field&>(static_cast<const Nodes*>(this)->field(name));
}

void Nodes::remove_field(const std::string& name)
{
  if( fields_.find(name)==fields_.end() ) {
    std::stringstream msg;
    msg << "Trying to remove field '"<<name<<"' from Nodes, but it is not present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_.erase(name);
}

void Nodes::resize( size_t size )
{
  FunctionSpace::resize(make_shape(size, shape(1)));

//  dof_ = size;
//  for( FieldMap::iterator it = fields_.begin(); it != fields_.end(); ++it )
//  {
//    Field& field = *it->second;
//    ArrayShape shape = field.shape();
//    shape[0] = size_;
//    field.resize(shape);
//  }
}

const Field& Nodes::field(size_t idx) const
{
  size_t c;
  for( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it )
  {
    if( idx == c )
    {
      const Field& field = *it->second;
      return field;
    }
    c++;
  }
}
Field& Nodes::field(size_t idx)
{
  return const_cast<Field&>(static_cast<const Nodes*>(this)->field(idx));
}


//-----------------------------------------------------------------------------


}  // namespace atlas

