/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

// file deepcode ignore CppMemoryLeak: static pointers for global registry are OK and will be cleaned up at end

#include <map>
#include <string>

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/output/Output.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

using atlas::FieldSet;
using atlas::FunctionSpace;
using atlas::Mesh;
using atlas::field::FieldImpl;
using atlas::field::FieldSetImpl;

namespace atlas {
namespace output {

static eckit::Mutex* local_mutex                        = nullptr;
static std::map<std::string, detail::OutputFactory*>* m = nullptr;
static pthread_once_t once                              = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    m           = new std::map<std::string, detail::OutputFactory*>();
}

namespace detail {

OutputImpl::OutputImpl() = default;

OutputImpl::~OutputImpl() = default;

}  // namespace detail

Output::Output(const std::string& key, std::ostream& stream, const eckit::Parametrisation& params):
    Handle(detail::OutputFactory::build(key, stream, params)) {}

/// Write mesh file
const Output& Output::write(const Mesh& m, const eckit::Parametrisation& c) const {
    get()->write(m, c);
    return *this;
}

/// Write field to file
const Output& Output::write(const Field& f, const eckit::Parametrisation& c) const {
    get()->write(f, c);
    return *this;
}

/// Write fieldset to file using FunctionSpace
const Output& Output::write(const FieldSet& f, const eckit::Parametrisation& c) const {
    get()->write(f, c);
    return *this;
}

/// Write field to file using Functionspace
const Output& Output::write(const Field& f, const FunctionSpace& fs, const eckit::Parametrisation& c) const {
    get()->write(f, fs, c);
    return *this;
}

/// Write fieldset to file using FunctionSpace
const Output& Output::write(const FieldSet& f, const FunctionSpace& fs, const eckit::Parametrisation& c) const {
    get()->write(f, fs, c);
    return *this;
}

namespace detail {

OutputFactory::OutputFactory(const std::string& name): name_(name) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ATLAS_ASSERT(m);
    if (m->find(name) != m->end()) {
        throw_Exception("Duplicate OutputFactory entry " + name);
    }

    (*m)[name] = this;
}

OutputFactory::~OutputFactory() {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    ATLAS_ASSERT(m);
    m->erase(name_);
}

void OutputFactory::list(std::ostream& out) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ATLAS_ASSERT(m);
    const char* sep = "";
    for (std::map<std::string, OutputFactory*>::const_iterator j = m->begin(); j != m->end(); ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}

const OutputImpl* OutputFactory::build(const std::string& name, std::ostream& stream) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ATLAS_ASSERT(m);
    std::map<std::string, OutputFactory*>::const_iterator j = m->find(name);

    Log::debug() << "Looking for OutputFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        Log::error() << "No OutputFactory for [" << name << "]" << std::endl;
        Log::error() << "OutputFactories are:" << std::endl;
        for (j = m->begin(); j != m->end(); ++j) {
            Log::error() << "   " << (*j).first << std::endl;
        }
        throw_Exception(std::string("No OutputFactory called ") + name);
    }

    return (*j).second->make(stream);
}

const OutputImpl* OutputFactory::build(const std::string& name, std::ostream& stream,
                                       const eckit::Parametrisation& param) {
    pthread_once(&once, init);
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ATLAS_ASSERT(m);
    std::map<std::string, OutputFactory*>::const_iterator j = m->find(name);

    Log::debug() << "Looking for OutputFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        Log::error() << "No OutputFactory for [" << name << "]" << std::endl;
        Log::error() << "OutputFactories are:" << std::endl;
        for (j = m->begin(); j != m->end(); ++j) {
            Log::error() << "   " << (*j).first << std::endl;
        }
        throw_Exception(std::string("No OutputFactory called ") + name);
    }

    return (*j).second->make(stream, param);
}

extern "C" {

void atlas__Output__delete(OutputImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_Output");
    delete This;
}

const OutputImpl* atlas__Output__create(const char* factory_key, std::ostream* stream,
                                        const eckit::Parametrisation* config) {
    ATLAS_ASSERT(config != nullptr, "Cannot access uninitialisd atlas_Config");
    const OutputImpl* output(nullptr);
    {
        Output o(std::string{factory_key}, *stream, *config);
        output = o.get();
        output->attach();
    }
    output->detach();
    return output;
}

void atlas__Output__write_mesh(const OutputImpl* This, Mesh::Implementation* mesh,
                               const eckit::Parametrisation* params) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_Output");
    Mesh m(mesh);
    This->write(m, *params);
}
void atlas__Output__write_fieldset(const OutputImpl* This, const FieldSetImpl* fieldset,
                                   const eckit::Parametrisation* config) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_Output");
    ATLAS_ASSERT(fieldset != nullptr, "Cannot access uninitialisd atlas_FieldSet");
    This->write(fieldset, *config);
}
void atlas__Output__write_field(const OutputImpl* This, const FieldImpl* field, const eckit::Parametrisation* config) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_Output");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialisd atlas_Field");
    ATLAS_ASSERT(config != nullptr, "Cannot access uninitialisd atlas_Config");
    This->write(field, *config);
}
void atlas__Output__write_fieldset_fs(const OutputImpl* This, const FieldSetImpl* fieldset,
                                      const functionspace::FunctionSpaceImpl* functionspace,
                                      const eckit::Parametrisation* params) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_Output");
    ATLAS_ASSERT(fieldset != nullptr, "Cannot access uninitialisd atlas_FieldSet");
    ATLAS_ASSERT(functionspace != nullptr, "Cannot access uninitialisd atlas_FunctionSpace");

    This->write(fieldset, functionspace, *params);
}
void atlas__Output__write_field_fs(const OutputImpl* This, const FieldImpl* field,
                                   const functionspace::FunctionSpaceImpl* functionspace,
                                   const eckit::Parametrisation* config) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_Output");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialisd atlas_Field");
    ATLAS_ASSERT(functionspace != nullptr, "Cannot access uninitialisd atlas_FunctionSpace");
    ATLAS_ASSERT(config != nullptr, "Cannot access uninitialisd atlas_Config");
    This->write(field, functionspace, *config);
}
}

}  // namespace detail
}  // namespace output
}  // namespace atlas
