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
#include "atlas/mesh/Nodes.h"
#include "atlas/Field.h"
#include "atlas/Parameters.h"

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

Nodes::Nodes(): size_(0)
{
  global_index_ = &add( Field::create<gidx_t>("glb_idx",   make_shape(size(),1)) );
  remote_index_ = &add( Field::create<int   >("remote_idx",make_shape(size(),1)) );
  partition_    = &add( Field::create<int   >("partition", make_shape(size(),1)) );
  lonlat_       = &add( Field::create<double>("lonlat",    make_shape(size(),2)) );
  ghost_        = &add( Field::create<int   >("ghost",     make_shape(size(),1)) );

  edge_connectivity_ = &add( "edge", new Connectivity() );
  cell_connectivity_ = &add( "cell", new Connectivity() );


  add( Field::create<int>("flags", make_shape(size(),1)) );


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
}

Nodes::Connectivity& Nodes::add( const std::string& name, Connectivity *connectivity )
{
  connectivities_[name] = connectivity;
  return *connectivity;
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
  return *field;
}

void Nodes::remove_field(const std::string& name)
{
  if( ! has_field(name) )
  {
    std::stringstream msg;
    msg << "Trying to remove field `"<<name<<"' in Nodes, but no field with this name is present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_.erase(name);
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

void Nodes::resize( size_t size )
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

const Field& Nodes::field(size_t idx) const
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
Field& Nodes::field(size_t idx)
{
  return const_cast<Field&>(static_cast<const Nodes*>(this)->field(idx));
}

void Nodes::print(std::ostream& os) const
{
    os << "Nodes[\n";
    os << "\t size=" << size() << ",\n";
    os << "\t fields=\n";
    for(size_t i = 0; i < nb_fields(); ++i)
    {
        os << "\t\t" << field(i);
        if( i != nb_fields()-1 )
          os << ",";
        os << "\n";
    }
    os << "]";
}

//-----------------------------------------------------------------------------

extern "C" {
int atlas__mesh__Nodes__size (Nodes* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->size();
  );
  return 0;
}
void atlas__mesh__Nodes__resize (Nodes* This, int size)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    This->resize(size);
  );
}
int atlas__mesh__Nodes__nb_fields (Nodes* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->nb_fields();
  );
  return 0;
}

void atlas__mesh__Nodes__add (Nodes* This, Field* field)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(field);
    This->add(field);
  );
}

void atlas__mesh__Nodes__remove_field (Nodes* This, char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    This->remove_field(std::string(name));
  );
}

int atlas__mesh__Nodes__has_field (Nodes* This, char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->has_field(std::string(name));
  );
  return 0;
}

Field* atlas__mesh__Nodes__field_by_name (Nodes* This, char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return &This->field(std::string(name));
  );
  return 0;
}

Field* atlas__mesh__Nodes__field_by_idx (Nodes* This, int idx)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return &This->field(idx);
  );
  return 0;
}

Metadata* atlas__mesh__Nodes__metadata(Nodes* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return &This->metadata();
  );
  return 0;
}

void atlas__mesh__Nodes__str (Nodes* This, char* &str, int &size)
{
  ATLAS_ERROR_HANDLING(
    std::stringstream ss;
    ss << *This;
    std::string s = ss.str();
    size = s.size();
    str = new char[size+1];
    strcpy(str,s.c_str());
  );
}

}

}  // namespace mesh
}  // namespace atlas

