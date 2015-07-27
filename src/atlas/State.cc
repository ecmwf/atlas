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

void State::initialize( const std::string& generator, const eckit::Parametrisation& params )
{
  eckit::ScopedPtr<StateGenerator> state_generator ( StateGeneratorFactory::build(generator, params) );
  state_generator->generate( *this, params );
}

//------------------------------------------------------------------------------------------------------

State::State()
{
}

State::State( const std::string& generator, const eckit::Parametrisation& params )
{
  initialize(generator,params);
}

const Metadata& State::metadata() const
{
  return metadata_;
}

Metadata& State::metadata()
{
  return metadata_;
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


void State::remove_field(const std::string& name)
{
  if( fields_.find(name)==fields_.end() ) {
    std::stringstream msg;
    msg << "Trying to remove field '"<<name<<"' from State, but it is not present in State.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_.erase(name);
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

void atlas__State__initialize(State* This, const char* generator, const eckit::Parametrisation* params)
{
  ASSERT( This );
  ASSERT( params );
  ATLAS_ERROR_HANDLING( This->initialize(std::string(generator),*params) );
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

Metadata* atlas__State__metadata (State* This)
{
  ASSERT( This );
  return &This->metadata();
}


}
//-----------------------------------------------------------------------------


}  // namespace atlas

