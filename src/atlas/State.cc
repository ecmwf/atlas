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
#include <string>
#include <map>

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/Field.h"
#include "atlas/Grid.h"
#include "atlas/Mesh.h"
#include "atlas/State.h"

namespace atlas {

namespace {

    static eckit::Mutex *local_mutex = 0;
    static std::map<std::string, StateFactory *> *m = 0;
    static pthread_once_t once = PTHREAD_ONCE_INIT;

    static void init() {
        local_mutex = new eckit::Mutex();
        m = new std::map<std::string, StateFactory *>();
    }

    template<typename T> void load_builder() { StateBuilder<T>("tmp"); }

    struct force_link {
        force_link()
        {
            // load_builder< A DERIVED TYPE >();
            // ...
        }
    };
}

State* State::create( const std::string& state_type, const eckit::Parametrisation& params)
{
  return StateFactory::build(state_type,params);
}

//------------------------------------------------------------------------------------------------------

State::State()
{
}

Field& State::add( Field* field )
{
  ASSERT( field != NULL );

  if( field->name().empty() )
  {
    std::stringstream new_name;
    new_name << "field_" << std::setw(5) << std::setfill('0') << fields_.size();
    ASSERT( !has_field(new_name.str() ) );
    field->name_ = new_name.str();
  }

  if( has_field(field->name()) ) {
    std::stringstream msg;
    msg << "Trying to add field '"<<field->name()<<"' to State, but State already has a field with this name.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_[field->name()] = eckit::SharedPtr<Field>(field);
  return *field;
}

Mesh& State::add( Mesh* mesh )
{
  ASSERT( mesh != NULL );
  ASSERT( meshes_.size() == 0 ); // Multiple meshes per state not yet supported
  meshes_[""] = eckit::SharedPtr<Mesh>(mesh);
  return *mesh;
}

Grid& State::add( Grid* grid )
{
  ASSERT( grid != NULL );
  ASSERT( grids_.size() == 0 ); // Multiple grids per state not yet supported
  grids_[""] = eckit::SharedPtr<Grid>(grid);
  return *grid;
}

const Field& State::field(const std::string& name) const
{
  if( ! has_field(name) )
  {
    std::stringstream msg;
    msg << "Trying to access field `"<<name<<"' in State, but no field with this name is present in State.";
    throw eckit::Exception(msg.str(),Here());
  }
  return *fields_.find(name)->second;
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
  FieldMap::const_iterator it = fields_.begin();
  for(size_t i = 0; i < idx; ++i) ++it;
  return *it->second;
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

  for( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it )
  {
    ret.push_back( it->first );
  }
  return ret;
}


const Mesh& State::mesh(const size_t idx) const
{
  ASSERT( idx < meshes_.size() );
  MeshMap::const_iterator it = meshes_.begin();
  for(size_t i = 0; i < idx; ++i) ++it;
  return *it->second;
}

Mesh& State::mesh(const size_t idx)
{
  ASSERT( idx < meshes_.size() );
  MeshMap::const_iterator it = meshes_.begin();
  for(size_t i = 0; i < idx; ++i) ++it;
  return *it->second;
}

const Grid& State::grid(const size_t idx) const
{
  ASSERT( idx < grids_.size() );
  GridMap::const_iterator it = grids_.begin();
  for(size_t i = 0; i < idx; ++i) ++it;
  return *it->second;
}

Grid& State::grid(const size_t idx)
{
  ASSERT( idx < grids_.size() );
  GridMap::const_iterator it = grids_.begin();
  for(size_t i = 0; i < idx; ++i) ++it;
  return *it->second;
}

const Mesh& State::mesh(const std::string& name) const { ASSERT( has_mesh(name) ); return *meshes_.find(name)->second; }
      Mesh& State::mesh(const std::string& name)       { ASSERT( has_mesh(name) ); return *meshes_[name]; }

const Grid& State::grid(const std::string& name) const { ASSERT( has_grid(name) ); return *grids_.find(name)->second; }
      Grid& State::grid(const std::string& name)       { ASSERT( has_grid(name) ); return *grids_[name]; }


void State::remove_field(const std::string& name)
{
  if( fields_.find(name)==fields_.end() ) {
    std::stringstream msg;
    msg << "Trying to remove field '"<<name<<"' from State, but it is not present in State.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_.erase(name);
}

void State::remove_mesh(const std::string& name)
{
  if( meshes_.find(name)==meshes_.end() ) {
    std::stringstream msg;
    msg << "Trying to remove mesh '"<<name<<"' from State, but it is not present in State.";
    throw eckit::Exception(msg.str(),Here());
  }
  meshes_.erase(name);
}

void State::remove_grid(const std::string& name)
{
  if( grids_.find(name)==grids_.end() ) {
    std::stringstream msg;
    msg << "Trying to remove grid '"<<name<<"' from State, but it is not present in State.";
    throw eckit::Exception(msg.str(),Here());
  }
  grids_.erase(name);
}

//-----------------------------------------------------------------------------


State* StateFactory::build(const std::string &name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, StateFactory *>::const_iterator j = m->find(name);

    eckit::Log::debug() << "Looking for StateFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        eckit::Log::error() << "No StateFactory for [" << name << "]" << std::endl;
        eckit::Log::error() << "StateFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j)
            eckit::Log::error() << "   " << (*j).first << std::endl;
        throw eckit::SeriousBug(std::string("No StateFactory called ") + name);
    }

    return (*j).second->make();
}

State* StateFactory::build(const std::string& name, const eckit::Parametrisation& param) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, StateFactory *>::const_iterator j = m->find(name);

    eckit::Log::debug() << "Looking for StateFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        eckit::Log::error() << "No StateFactory for [" << name << "]" << std::endl;
        eckit::Log::error() << "StateFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j)
            eckit::Log::error() << "   " << (*j).first << std::endl;
        throw eckit::SeriousBug(std::string("No StateFactory called ") + name);
    }

    return (*j).second->make(param);
}




//-----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"{

}
//-----------------------------------------------------------------------------


}  // namespace atlas

