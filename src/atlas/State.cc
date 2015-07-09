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
    static std::map<std::string, StateGeneratorFactory *> *m = 0;
    static pthread_once_t once = PTHREAD_ONCE_INIT;

    static void init() {
        local_mutex = new eckit::Mutex();
        m = new std::map<std::string, StateGeneratorFactory *>();
    }

    template<typename T> void load_builder() { StateGeneratorBuilder<T>("tmp"); }

    struct force_link {
        force_link()
        {
            // load_builder< A DERIVED TYPE >();
            // ...
        }
    };
}

State* State::create( const std::string& _state_generator, const eckit::Parametrisation& params)
{
  eckit::ScopedPtr<StateGenerator> state_generator ( StateGeneratorFactory::build(_state_generator, params) );
  State* state = new State;
  state_generator->generate( *state, params );
  return state;
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
  if( meshes_.size() != 0 )
    throw eckit::NotImplemented("Multiple meshes per state not yet supported",Here());
  meshes_[""] = eckit::SharedPtr<Mesh>(mesh);
  return *mesh;
}

Grid& State::add( Grid* grid )
{
  ASSERT( grid != NULL );
  if( grids_.size() != 0 )
    throw eckit::NotImplemented("Multiple grids per state not yet supported",Here());
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

StateGenerator::StateGenerator( const eckit::Parametrisation& )
{}

StateGenerator::~StateGenerator()
{}

StateGenerator* StateGeneratorFactory::build(const std::string& name, const eckit::Parametrisation& param) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, StateGeneratorFactory *>::const_iterator j = m->find(name);

    eckit::Log::debug() << "Looking for StateGeneratorFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        eckit::Log::error() << "No StateGeneratorFactory for [" << name << "]" << std::endl;
        eckit::Log::error() << "StateFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j)
            eckit::Log::error() << "   " << (*j).first << std::endl;
        throw eckit::SeriousBug(std::string("No StateGeneratorFactory called ") + name);
    }

    return (*j).second->make(param);
}

void StateGeneratorFactory::list(std::ostream& out) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    const char* sep = "";
    for (std::map<std::string, StateGeneratorFactory *>::const_iterator j = m->begin() ; j != m->end() ; ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}

bool StateGeneratorFactory::has(const std::string& name)
{
  pthread_once(&once, init);

  eckit::AutoLock<eckit::Mutex> lock(local_mutex);

  static force_link static_linking;

  return ( m->find(name) != m->end() );
}

StateGeneratorFactory::StateGeneratorFactory(const std::string &name):
    name_(name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}


StateGeneratorFactory::~StateGeneratorFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}



//-----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C"{

State* atlas__State__new()
{
  return new State;
}

State* atlas__State__create(const char* factory, const eckit::Parametrisation* params)
{
  ASSERT( params );
  State* state(0);
  ATLAS_ERROR_HANDLING( state = State::create(std::string(factory),*params) );
  return state;
}

void atlas__State__delete (State* This)
{
  ASSERT( This );
  delete This;
}

void atlas__State__add_field (State* This, Field* field)
{
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->add(field); );
}

void atlas__State__remove_field (State* This, const char* name)
{
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->remove_field(name); );
}

int atlas__State__has_field (State* This, const char* name)
{
  ASSERT( This );
  int has_field(0);
  ATLAS_ERROR_HANDLING( has_field = This->has_field(name); );
  return has_field;
}

Field* atlas__State__field_by_name (State* This, const char* name)
{
  ASSERT( This );
  Field* field(0);
  ATLAS_ERROR_HANDLING ( field = &This->field( std::string(name) ); );
  return field;
}

Field* atlas__State__field_by_index (State* This, int index)
{
  ASSERT( This );
  Field* field(0);
  ATLAS_ERROR_HANDLING( field = &This->field( index ) );
  return field;
}

int atlas__State__nb_fields(const State* This)
{
  ASSERT( This );
  int nb_fields(0);
  ATLAS_ERROR_HANDLING( nb_fields = This->nb_fields(); );
  return nb_fields;
}

void atlas__State__add_grid (State* This, Grid* grid)
{
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->add(grid); );
}

void atlas__State__remove_grid (State* This, const char* name)
{
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->remove_grid(name); );
}

int atlas__State__has_grid (State* This, const char* name)
{
  ASSERT( This );
  int has_grid(0);
  ATLAS_ERROR_HANDLING( has_grid = This->has_grid(name); );
  return has_grid;
}

Grid* atlas__State__grid_by_name (State* This, const char* name)
{
  ASSERT( This );
  Grid* grid(0);
  ATLAS_ERROR_HANDLING( grid = &This->grid( std::string(name) ); );
  return grid;
}

Grid* atlas__State__grid_by_index (State* This, int index)
{
  ASSERT( This );
  Grid* grid(0);
  ATLAS_ERROR_HANDLING( grid = &This->grid( index ); );
  return grid;
}

int atlas__State__nb_grids(const State* This)
{
  ASSERT( This );
  int nb_grids(0);
  ATLAS_ERROR_HANDLING( nb_grids = This->nb_grids(); );
  return nb_grids;
}

void atlas__State__add_mesh (State* This, Mesh* mesh)
{
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->add(mesh); );
}

void atlas__State__remove_mesh (State* This, const char* name)
{
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->remove_mesh(name); );
}

int atlas__State__has_mesh (State* This, const char* name)
{
  ASSERT( This );
  int has_mesh(0);
  ATLAS_ERROR_HANDLING( has_mesh = This->has_mesh(name); );
  return has_mesh;
}

Mesh* atlas__State__mesh_by_name (State* This, const char* name)
{
  ASSERT( This );
  Mesh* mesh(0);
  ATLAS_ERROR_HANDLING( mesh = &This->mesh( std::string(name) ); );
  return mesh;
}

Mesh* atlas__State__mesh_by_index (State* This, int index)
{
  ASSERT( This );
  Mesh* mesh(0);
  ATLAS_ERROR_HANDLING( mesh = &This->mesh( index ); );
  return mesh;
}

int atlas__State__nb_meshes(const State* This)
{
  int nb_meshes(0);
  ATLAS_ERROR_HANDLING( nb_meshes = This->nb_meshes(); );
  return nb_meshes;
}

}
//-----------------------------------------------------------------------------


}  // namespace atlas

