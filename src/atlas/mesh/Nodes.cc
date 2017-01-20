/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/array/MakeView.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/field/Field.h"
#include "atlas/internals/Parameters.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

Nodes::Nodes(): size_(0)
{
  global_index_ = &add( field::Field::create<gidx_t>("glb_idx",   array::make_shape(size())) );
  remote_index_ = &add( field::Field::create<int   >("remote_idx",array::make_shape(size())) );
  partition_    = &add( field::Field::create<int   >("partition", array::make_shape(size())) );
  lonlat_       = &add( field::Field::create<double>("lonlat",    array::make_shape(size(),2)) );
  ghost_        = &add( field::Field::create<int   >("ghost",     array::make_shape(size())) );

  edge_connectivity_ = &add( new Connectivity("edge") );
  cell_connectivity_ = &add( new Connectivity("cell") );


  add( field::Field::create<int>("flags", array::make_shape(size())) );


  array::ArrayView<gidx_t,1> glb_idx = array::make_view<gidx_t,1>( global_index() );
  array::ArrayView<int   ,1> part = array::make_view<int,1>( partition() );
  array::ArrayView<int   ,1> flags = array::make_view<int,1>( field("flags") );

  for(size_t n=0; n<size(); ++n)
  {
    glb_idx(n) = 1+n;
    part(n) = parallel::mpi::comm().rank();
    flags(n) = 0;
  }
}

Nodes::Connectivity& Nodes::add( Connectivity *connectivity )
{
  connectivities_[connectivity->name()] = connectivity;
  return *connectivity;
}

field::Field& Nodes::add( field::Field* field )
{
  ASSERT( field != NULL );
  ASSERT( ! field->name().empty() );

  if( has_field(field->name()) ) {
    std::stringstream msg;
    msg << "Trying to add field '"<<field->name()<<"' to Nodes, but Nodes already has a field with this name.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_[field->name()] = eckit::SharedPtr<field::Field>(field);
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


const field::Field& Nodes::field(const std::string& name) const
{
  if( ! has_field(name) )
  {
    std::stringstream msg;
    msg << "Trying to access field `"<<name<<"' in Nodes, but no field with this name is present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  return *fields_.find(name)->second;
}

field::Field& Nodes::field(const std::string& name)
{
  return const_cast<field::Field&>(static_cast<const Nodes*>(this)->field(name));
}

void Nodes::resize( size_t size )
{
  size_t previous_size = size_;
  size_ = size;
  for( FieldMap::iterator it = fields_.begin(); it != fields_.end(); ++it )
  {
    field::Field& field = *it->second;
    array::ArrayShape shape = field.shape();
    shape[0] = size_;
    field.resize(shape);
  }

  array::ArrayView<gidx_t,1> glb_idx = array::make_view<gidx_t,1>( global_index() );
  array::ArrayView<int   ,1> part    = array::make_view<int,   1>( partition() );
  array::ArrayView<int   ,1> flags   = array::make_view<int,   1>( field("flags") );

  for(size_t n=previous_size; n<size_; ++n)
  {
    glb_idx(n) = 1+n;
    part(n) = parallel::mpi::comm().rank();
    flags(n) = 0;
  }
}

const field::Field& Nodes::field(size_t idx) const
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
  eckit::SeriousBug("Should not be here!",Here());
  static field::Field* ret;
  return *ret;
}
field::Field& Nodes::field(size_t idx)
{
  return const_cast<field::Field&>(static_cast<const Nodes*>(this)->field(idx));
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


size_t Nodes::footprint() const {
  size_t size = sizeof(*this);
  for( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it ) {
    size += (*it).second->footprint();
  }
  for( ConnectivityMap::const_iterator it = connectivities_.begin(); it != connectivities_.end(); ++it ) {
    size += (*it).second->footprint();
  }
  size += metadata_.footprint();
  return size;
}


const IrregularConnectivity& Nodes::connectivity(const std::string& name) const
{
  if( (connectivities_.find(name) == connectivities_.end()) )
  {
    std::stringstream msg;
    msg << "Trying to access connectivity `"<<name<<"' in Nodes, but no connectivity with this name is present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  return *connectivities_.find(name)->second;
}
IrregularConnectivity& Nodes::connectivity(const std::string& name)
{
  if( (connectivities_.find(name) == connectivities_.end()) )
  {
    std::stringstream msg;
    msg << "Trying to access connectivity `"<<name<<"' in Nodes, but no connectivity with this name is present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  return *connectivities_.find(name)->second;
}

void Nodes::cloneToDevice() const {
  std::for_each(fields_.begin(), fields_.end(), [](const FieldMap::value_type& v){ v.second->cloneToDevice();});
  std::for_each(connectivities_.begin(), connectivities_.end(), [](const ConnectivityMap::value_type& v){ v.second->cloneToDevice();});
}

void Nodes::cloneFromDevice() const {
  std::for_each(fields_.begin(), fields_.end(), [](const FieldMap::value_type& v){ v.second->cloneFromDevice();});
  std::for_each(connectivities_.begin(), connectivities_.end(), [](const ConnectivityMap::value_type& v){ v.second->cloneFromDevice();});
}

void Nodes::syncHostDevice() const {
  std::for_each(fields_.begin(), fields_.end(), [](const FieldMap::value_type& v){ v.second->syncHostDevice();});
  std::for_each(connectivities_.begin(), connectivities_.end(), [](const ConnectivityMap::value_type& v){ v.second->syncHostDevice();});
}

//-----------------------------------------------------------------------------

extern "C" {

Nodes* atlas__mesh__Nodes__create()
{
  Nodes* nodes(0);
  ATLAS_ERROR_HANDLING( nodes = new Nodes() );
  return nodes;
}

void atlas__mesh__Nodes__delete (Nodes* This)
{
  ATLAS_ERROR_HANDLING( delete This );
}


size_t atlas__mesh__Nodes__size (Nodes* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->size();
  );
  return 0;
}
void atlas__mesh__Nodes__resize (Nodes* This, size_t size)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    This->resize(size);
  );
}
size_t atlas__mesh__Nodes__nb_fields (Nodes* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->nb_fields();
  );
  return 0;
}

void atlas__mesh__Nodes__add_field (Nodes* This, field::Field* field)
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

field::Field* atlas__mesh__Nodes__field_by_name (Nodes* This, char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return &This->field(std::string(name));
  );
  return 0;
}

field::Field* atlas__mesh__Nodes__field_by_idx (Nodes* This, size_t idx)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return &This->field(idx);
  );
  return 0;
}

util::Metadata* atlas__mesh__Nodes__metadata(Nodes* This)
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

IrregularConnectivity* atlas__mesh__Nodes__edge_connectivity(Nodes* This)
{
  IrregularConnectivity* connectivity(0);
  ATLAS_ERROR_HANDLING( connectivity = &This->edge_connectivity() );
  return connectivity;
}

IrregularConnectivity* atlas__mesh__Nodes__cell_connectivity(Nodes* This)
{
  IrregularConnectivity* connectivity(0);
  ATLAS_ERROR_HANDLING( connectivity = &This->cell_connectivity() );
  return connectivity;
}

IrregularConnectivity* atlas__mesh__Nodes__connectivity (Nodes* This, char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return &This->connectivity(std::string(name));
  );
  return 0;
}


void atlas__mesh__Nodes__add_connectivity (Nodes* This, IrregularConnectivity* connectivity)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(connectivity);
    This->add(connectivity);
  );
}

field::Field* atlas__mesh__Nodes__lonlat(Nodes* This)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->lonlat();
  );
  return field;
}

field::Field* atlas__mesh__Nodes__global_index(Nodes* This)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->global_index();
  );
  return field;
}

field::Field* atlas__mesh__Nodes__remote_index(Nodes* This)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->remote_index();
  );
  return field;
}

field::Field* atlas__mesh__Nodes__partition(Nodes* This)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->partition();
  );
  return field;
}

field::Field* atlas__mesh__Nodes__ghost(Nodes* This)
{
  field::Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->ghost();
  );
  return field;
}
}

}  // namespace mesh
}  // namespace atlas

