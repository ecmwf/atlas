/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */
#include <sstream>
#include <map>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/log/Log.h"

#include "atlas/field/FieldCreator.h"
#include "atlas/Field.h"
#include "atlas/Grid.h"
#include "atlas/FunctionSpace.h"

namespace {

    static eckit::Mutex *local_mutex = 0;
    static std::map<std::string, atlas::field::FieldCreatorFactory *> *m = 0;
    static pthread_once_t once = PTHREAD_ONCE_INIT;

    static void init() {
        local_mutex = new eckit::Mutex();
        m = new std::map<std::string, atlas::field::FieldCreatorFactory *>();
    }
}

namespace atlas {
namespace field {

// ------------------------------------------------------------------

FieldCreator::FieldCreator()
{
}

FieldCreator::~FieldCreator()
{
}

FieldCreatorFactory::FieldCreatorFactory(const std::string &name):
    name_(name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}


FieldCreatorFactory::~FieldCreatorFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}


void FieldCreatorFactory::list(std::ostream& out) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    const char* sep = "";
    for (std::map<std::string, FieldCreatorFactory *>::const_iterator j = m->begin() ; j != m->end() ; ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}


FieldCreator *FieldCreatorFactory::build(const std::string &name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    std::map<std::string, FieldCreatorFactory *>::const_iterator j = m->find(name);

    eckit::Log::debug() << "Looking for FieldCreatorFactory [" << name << "]" << '\n';

    if (j == m->end()) {
        eckit::Log::error() << "No FieldCreatorFactory for [" << name << "]" << '\n';
        eckit::Log::error() << "FieldCreatorFactories are:" << '\n';
        for (j = m->begin() ; j != m->end() ; ++j)
            eckit::Log::error() << "   " << (*j).first << '\n';
        throw eckit::SeriousBug(std::string("No FieldCreatorFactory called ") + name);
    }

    return (*j).second->make();
}

FieldCreator *FieldCreatorFactory::build(const std::string& name, const eckit::Parametrisation& param) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    std::map<std::string, FieldCreatorFactory *>::const_iterator j = m->find(name);

    eckit::Log::debug() << "Looking for FieldCreatorFactory [" << name << "]" << '\n';

    if (j == m->end()) {
        eckit::Log::error() << "No FieldCreatorFactory for [" << name << "]" << '\n';
        eckit::Log::error() << "FieldCreatorFactories are:" << '\n';
        for (j = m->begin() ; j != m->end() ; ++j)
            eckit::Log::error() << "   " << (*j).first << '\n';
        throw eckit::SeriousBug(std::string("No FieldCreatorFactory called ") + name);
    }

    return (*j).second->make(param);
}

} // namespace field
} // namespace atlas

