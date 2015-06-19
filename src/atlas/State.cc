/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iomanip>
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/Field.h"
#include "atlas/Grid.h"
#include "atlas/Mesh.h"
#include "atlas/State.h"

namespace atlas {

//------------------------------------------------------------------------------------------------------

State::State()
{
}

State::State( const eckit::Parametrisation& )
{
}

Field& State::add( Field* field )
{
  ASSERT( field != NULL );

  if( field->name().empty() )
  {
    std::stringstream new_name;
    new_name << "field_" << std::setw(5) << std::setfill('0') << fields_.size();
    ASSERT( ! fields_.has(new_name.str()) );
    field->name_ = new_name.str();
  }

  if( fields_.has(field->name()) ) {
    std::stringstream msg;
    msg << "Trying to add field '"<<field->name()<<"' to State, but State already has a field with this name.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_.insert(field->name(),eckit::SharedPtr<Field>(field));
  fields_.sort();
  return *field;
}

void State::set_name( Field& field )
{
}


Mesh& State::add( Mesh* mesh )
{
  ASSERT( mesh != NULL );
  meshes_.insert("",eckit::SharedPtr<Mesh>(mesh));
  meshes_.sort();
  return *mesh;
}

Grid& State::add( Grid* grid )
{
  ASSERT( grid != NULL );
  grids_.insert("",eckit::SharedPtr<Grid>(grid));
  grids_.sort();
  return *grid;
}

const Field& State::field(const std::string& name) const
{
  if( !fields_.has(name) )
  {
    std::stringstream msg;
    msg << "Trying to access field `"<<name<<"' in State, but no field with this name is present in State.";
    throw eckit::Exception(msg.str(),Here());
  }
  return *fields_[name];
}

Field& State::field(const std::string& name)
{
  return const_cast<Field&>(static_cast<const State*>(this)->field(name));
}

const Field& State::field(const size_t idx) const
{
  if( idx >= fields_.size() )
  {
    std::stringstream msg;
    msg << "Trying to access field in State with index "<<idx<<", but there exist only "<<fields_.size()<<" fields in State.";
    throw eckit::Exception(msg.str(),Here());
  }
  return *fields_[idx];
}

Field& State::field(const size_t idx)
{
  return const_cast<Field&>(static_cast<const State*>(this)->field(idx));
}

std::vector< std::string > State::field_names() const
{
  std::vector< std::string > ret;
  if (fields_.size())
    ret.reserve(fields_.size());

  for( FieldMap::const_iterator it = fields_.cbegin(); it != fields_.cend(); ++it )
  {
    ret.push_back( it->key );
  }
  return ret;
}


const Mesh& State::mesh(const size_t idx) const { ASSERT( idx < meshes_.size() ); return *meshes_[idx]; }
      Mesh& State::mesh(const size_t idx)       { ASSERT( idx < meshes_.size() ); return *meshes_[idx]; }

const Grid& State::grid(const size_t idx) const { ASSERT( idx < grids_.size() ); return *grids_[idx]; }
      Grid& State::grid(const size_t idx)       { ASSERT( idx < grids_.size() ); return *grids_[idx]; }

const Mesh& State::mesh(const std::string& name) const { ASSERT( meshes_.has(name) ); return *meshes_[name]; }
      Mesh& State::mesh(const std::string& name)       { ASSERT( meshes_.has(name) ); return *meshes_[name]; }

const Grid& State::grid(const std::string& name) const { ASSERT( grids_.has(name) ); return *grids_[name]; }
      Grid& State::grid(const std::string& name)       { ASSERT( grids_.has(name) ); return *grids_[name]; }

#if 0 // needs DenseMap to implment erase

void State::remove_field(const std::string& name)
{
  if( !fields_.has(name) ) {
    std::stringstream msg;
    msg << "Trying to remove field '"<<name<<"' from State, but it is not present in State.";
    throw eckit::Exception(msg.str(),Here());
  }
  eckit::SharedPtr<Field> field = fields_[name];
  fields_.erase(name);
  fields_.sort();
}

void State::remove_mesh(const std::string& name)
{
  if( !meshes_.has(name) ) {
    std::stringstream msg;
    msg << "Trying to remove mesh "<<name<<" from State, but it is not present in State.";
    throw eckit::Exception(msg.str(),Here());
  }
  eckit::SharedPtr<Mesh> mesh = meshes_[name];
  meshes_.erase(name);
  meshes_.sort();
}

void State::remove_grid(const std::string& name)
{
  if( !grids_.has(name) ) {
    std::stringstream msg;
    msg << "Trying to remove grid "<<name<<" from State, but it is not present in State.";
    throw eckit::Exception(msg.str(),Here());
  }
  eckit::SharedPtr<Grid> grid = grids_[name];
  grids_.erase(name);
  grids_.sort();
}

#endif

//-----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"{

}
//-----------------------------------------------------------------------------


}  // namespace atlas

