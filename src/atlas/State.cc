/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

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
  if( !field->name().empty() && fields_.has(field->name()) ) {
    std::stringstream msg;
    msg << "Trying to add field '"<<field->name()<<"' to State, but State already has a field with this name.";
    throw eckit::Exception(msg.str(),Here());
  }
  if( field->name().empty() )
  {
    // The fact there is no name, means it is not important to reaccess it by name
    std::stringstream new_name;
    new_name << "field_" << fields_.size();
    ASSERT( ! fields_.has(new_name.str()) );
    fields_.insert(new_name.str(),eckit::SharedPtr<Field>(field));
  }
  else
  {
    fields_.insert(field->name(),eckit::SharedPtr<Field>(field));
  }
  fields_.sort();
  return *field;
}

Mesh& State::add( Mesh* mesh )
{
  ASSERT( mesh != NULL );
  meshes_.push_back( eckit::SharedPtr<Mesh>(mesh) );
  return *mesh;
}

Grid& State::add( Grid* grid )
{
  ASSERT( grid != NULL );
  grids_.push_back( eckit::SharedPtr<Grid>(grid) );
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

const Mesh& State::mesh(const size_t idx) const { ASSERT( idx < meshes_.size() ); return *meshes_[idx]; }
      Mesh& State::mesh(const size_t idx)       { ASSERT( idx < meshes_.size() ); return *meshes_[idx]; }

const Grid& State::grid(const size_t idx) const { ASSERT( idx < grids_.size() ); return *grids_[idx]; }
      Grid& State::grid(const size_t idx)       { ASSERT( idx < grids_.size() ); return *grids_[idx]; }


//-----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"{

}
//-----------------------------------------------------------------------------


}  // namespace atlas

