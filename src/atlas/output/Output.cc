/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <map>
#include <string>

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/output/Output.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"

using atlas::field::FieldImpl;
using atlas::FieldSet;
using atlas::field::FieldSetImpl;
using atlas::mesh::Mesh;
using atlas::functionspace::FunctionSpace;
using eckit::Parametrisation;

namespace atlas {
namespace output {

namespace {

    static eckit::Mutex *local_mutex = 0;
    static std::map<std::string, OutputFactory *> *m = 0;
    static pthread_once_t once = PTHREAD_ONCE_INIT;

    static void init() {
        local_mutex = new eckit::Mutex();
        m = new std::map<std::string, OutputFactory *>();
    }

    template<typename T> void load_builder() { OutputBuilder<T>("tmp"); }

    struct force_link {
        force_link()
        {
            load_builder<output::detail::Gmsh>();
        }
    };

}

OutputImpl::OutputImpl() {
}

OutputImpl::~OutputImpl() {
}

Output::Output() :
    output_( nullptr ) {
}

Output::Output( const output_t* output ) :
    output_( output ) {
}

Output::Output( const Output& output ) :
    output_( output.output_ ) {
}


Output::Output(const std::string &key, Stream& stream, const eckit::Parametrisation &params) :
    output_( OutputFactory::build(key,stream,params) ) {
}


/// Write mesh file
void Output::write(
    const mesh::Mesh& m,
    const eckit::Parametrisation& c ) const {
  return output_->write(m,c);
}

/// Write field to file
void Output::write(
    const Field& f,
    const eckit::Parametrisation& c ) const {
  return output_->write(f,c);
}

/// Write fieldset to file using FunctionSpace
void Output::write(
    const FieldSet& f,
    const eckit::Parametrisation& c ) const {
  return output_->write(f,c);
}

/// Write field to file using Functionspace
void Output::write(
    const Field& f,
    const functionspace::FunctionSpace& fs,
    const eckit::Parametrisation& c ) const {
  return output_->write(f,fs,c);
}

/// Write fieldset to file using FunctionSpace
void Output::write(
    const FieldSet& f,
    const functionspace::FunctionSpace& fs,
    const eckit::Parametrisation& c ) const {
  return output_->write(f,fs,c);
}



OutputFactory::OutputFactory(const std::string &name):
    name_(name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}


OutputFactory::~OutputFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}


void OutputFactory::list(std::ostream& out) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    const char* sep = "";
    for (std::map<std::string, OutputFactory *>::const_iterator j = m->begin() ; j != m->end() ; ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}


const OutputImpl *OutputFactory::build(const std::string &name, Stream& stream) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, OutputFactory *>::const_iterator j = m->find(name);

    Log::debug<Atlas>() << "Looking for OutputFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        Log::error() << "No OutputFactory for [" << name << "]" << std::endl;
        Log::error() << "OutputFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j)
            Log::error() << "   " << (*j).first << std::endl;
        throw eckit::SeriousBug(std::string("No OutputFactory called ") + name);
    }

    return (*j).second->make(stream);
}

const OutputImpl *OutputFactory::build(const std::string& name, Stream& stream, const eckit::Parametrisation& param) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, OutputFactory *>::const_iterator j = m->find(name);

    Log::debug<Atlas>() << "Looking for OutputFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        Log::error() << "No OutputFactory for [" << name << "]" << std::endl;
        Log::error() << "OutputFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j)
            Log::error() << "   " << (*j).first << std::endl;
        throw eckit::SeriousBug(std::string("No OutputFactory called ") + name);
    }

    return (*j).second->make(stream,param);
}

extern "C" {

void atlas__Output__delete(OutputImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

const OutputImpl* atlas__Output__create(const char* factory_key, Stream* stream, const eckit::Parametrisation* params)
{
  const OutputImpl* output(0);
  ATLAS_ERROR_HANDLING (
    // ASSERT(stream);
    ASSERT(params);
    {
       Output o( std::string{factory_key}, *stream, *params );
       output = o.get();
       output->attach();
    }
    output->detach();
  );
  return output;
}

void atlas__Output__write_mesh(const OutputImpl* This, Mesh::Implementation* mesh, const Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(mesh);
    ASSERT(params);
    mesh::Mesh m(mesh);
    This->write(m,*params);
  );
}
void atlas__Output__write_fieldset(const OutputImpl* This, const FieldSetImpl* fieldset, const Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(fieldset);
    ASSERT(params);
    This->write(fieldset,*params);
  );
}
void atlas__Output__write_field(const OutputImpl* This, const FieldImpl* field, const Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(field);
    ASSERT(params);
    This->write(field,*params);
  );
}
void atlas__Output__write_fieldset_fs(const OutputImpl* This, const FieldSetImpl* fieldset, const functionspace::FunctionSpaceImpl* functionspace, const Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(fieldset);
    ASSERT(functionspace);
    ASSERT(params);
    This->write(fieldset,functionspace,*params);
  );
}
void atlas__Output__write_field_fs(const OutputImpl* This, const FieldImpl* field, const functionspace::FunctionSpaceImpl* functionspace, const Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(field);
    ASSERT(functionspace);
    ASSERT(params);
    This->write(field,functionspace,*params);
  );
}

}

} // namespace output
} // namespace atlas

