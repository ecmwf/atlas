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

using atlas::field::Field;
using atlas::field::FieldSet;
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
            load_builder<output::Gmsh>();
        }
    };

}

Output* Output::create(const std::string &key, Stream& stream, const eckit::Parametrisation &params)
{
  return OutputFactory::build(key,stream,params);
}

Output::Output()
{
}

Output::~Output() {
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


Output *OutputFactory::build(const std::string &name, Stream& stream) {

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

Output *OutputFactory::build(const std::string& name, Stream& stream, const eckit::Parametrisation& param) {

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

void atlas__Output__delete(Output* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

Output* atlas__Output__create(const char* factory_key, Stream* stream, const eckit::Parametrisation* params)
{
  Output* Output(0);
  ATLAS_ERROR_HANDLING (
    // ASSERT(stream);
    ASSERT(params);
    Output = Output::create(std::string(factory_key),*stream,*params);
  );
  return Output;
}

void atlas__Output__write_mesh(const Output* This, const Mesh* mesh, const Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(mesh);
    ASSERT(params);
    This->write(*mesh,*params);
  );
}
void atlas__Output__write_fieldset(const Output* This, const FieldSet* fieldset, const Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(fieldset);
    ASSERT(params);
    This->write(*fieldset,*params);
  );
}
void atlas__Output__write_field(const Output* This, const Field* field, const Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(field);
    ASSERT(params);
    This->write(*field,*params);
  );
}
void atlas__Output__write_fieldset_fs(const Output* This, const FieldSet* fieldset, const FunctionSpace* functionspace, const Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(fieldset);
    ASSERT(functionspace);
    ASSERT(params);
    This->write(*fieldset,*functionspace,*params);
  );
}
void atlas__Output__write_field_fs(const Output* This, const Field* field, const FunctionSpace* functionspace, const Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(field);
    ASSERT(functionspace);
    ASSERT(params);
    This->write(*field,*functionspace,*params);
  );
}

}

} // namespace output
} // namespace atlas

