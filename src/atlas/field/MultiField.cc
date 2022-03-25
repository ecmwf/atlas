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

#include "MultiField.h"

#include <iomanip>
#include <map>
#include <memory>
#include <string>

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"

#include "atlas/field/Field.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace field {

namespace {

static eckit::Mutex* local_mutex                           = nullptr;
static std::map<std::string, MultiFieldCreatorFactory*>* m = nullptr;
static pthread_once_t once                                 = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    m           = new std::map<std::string, MultiFieldCreatorFactory*>();
}

template <typename T>
void load_builder() {
    MultiFieldCreatorBuilder<T>("tmp");
}

void force_link() {
    static struct Link { Link() = default; } link;
    [](const Link&) {}(link);  // disable unused warnings
}

}  // namespace

void MultiField::initialize(const std::string& generator, const eckit::Parametrisation& params) {
    std::unique_ptr<MultiFieldCreator> FieldPool_generator(MultiFieldCreatorFactory::build(generator, params));
    FieldPool_generator->generate(*this, params);
}

array::Array& MultiField::allocate(array::DataType datatype, array::ArraySpec&& spec) {
    array_.reset(array::Array::create(datatype, std::move(spec)));
    return *array_;
}

//------------------------------------------------------------------------------------------------------

MultiField::MultiField() = default;

MultiField::MultiField(const std::string& generator, const eckit::Parametrisation& params) {
    initialize(generator, params);
}

const util::Metadata& MultiField::metadata() const {
    return metadata_;
}

util::Metadata& MultiField::metadata() {
    return metadata_;
}

std::vector<std::string> MultiField::field_names() const {
    std::vector<std::string> ret;
    if (fields_.size()) {
        ret.reserve(fields_.size());
    }

    for (auto it = field_index_.begin(); it != field_index_.end(); ++it) {
        ret.push_back(it->first);
    }
    return ret;
}

//-----------------------------------------------------------------------------

MultiFieldCreator::MultiFieldCreator(const eckit::Parametrisation&) {}

MultiFieldCreator::~MultiFieldCreator() = default;

MultiFieldCreator* MultiFieldCreatorFactory::build(const std::string& name, const eckit::Parametrisation& param) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    force_link();

    std::map<std::string, MultiFieldCreatorFactory*>::const_iterator j = m->find(name);

    Log::debug() << "Looking for MultiFieldCreatorFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        Log::error() << "No MultiFieldCreatorFactory for [" << name << "]" << std::endl;
        Log::error() << "FieldPoolFactories are:" << std::endl;
        for (j = m->begin(); j != m->end(); ++j) {
            Log::error() << "   " << (*j).first << std::endl;
        }
        throw_Exception(std::string("No MultiFieldCreatorFactory called ") + name, Here());
    }

    return (*j).second->make(param);
}

void MultiFieldCreatorFactory::list(std::ostream& out) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    force_link();

    const char* sep = "";
    for (std::map<std::string, MultiFieldCreatorFactory*>::const_iterator j = m->begin(); j != m->end(); ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}

bool MultiFieldCreatorFactory::has(const std::string& name) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    force_link();

    return (m->find(name) != m->end());
}

MultiFieldCreatorFactory::MultiFieldCreatorFactory(const std::string& name): name_(name) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ATLAS_ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}

MultiFieldCreatorFactory::~MultiFieldCreatorFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}

//-----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {

MultiField* atlas__FieldPool__new() {
    return new MultiField;
}

void atlas__FieldPool__initialize(MultiField* This, const char* generator, const eckit::Parametrisation* params) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldPool");
    This->initialize(std::string(generator), *params);
}

void atlas__FieldPool__delete(MultiField* This) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldPool");
    delete This;
}

int atlas__FieldPool__has(MultiField* This, const char* name) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldPool");
    return This->has(name);
}

FieldImpl* atlas__FieldPool__field_by_name(MultiField* This, const char* name) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldPool");
    return This->field(std::string(name)).get();
}

FieldImpl* atlas__FieldPool__field_by_index(MultiField* This, idx_t index) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldPool");
    return This->field(index).get();
}

idx_t atlas__FieldPool__size(const MultiField* This) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldPool");
    return This->size();
}

util::Metadata* atlas__FieldPool__metadata(MultiField* This) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldPool");
    return &This->metadata();
}
}
//-----------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
