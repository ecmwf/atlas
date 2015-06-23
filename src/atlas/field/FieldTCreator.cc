/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "atlas/field/FieldTCreator.h"
#include "atlas/field/FieldT.h"


namespace {

    static eckit::Mutex *local_mutex = 0;
    static std::map<std::string, atlas::field::FieldTCreatorFactory *> *m = 0;
    static pthread_once_t once = PTHREAD_ONCE_INIT;

    static void init() {
        local_mutex = new eckit::Mutex();
        m = new std::map<std::string, atlas::field::FieldTCreatorFactory *>();
    }
}

namespace atlas {
namespace field {

namespace {

template<typename T> void load_builder() { FieldTCreatorBuilder<T>("tmp"); }

struct force_link {
    force_link()
    {
        load_builder< FieldTCreatorT<int> >();
        load_builder< FieldTCreatorT<long> >();
        load_builder< FieldTCreatorT<float> >();
        load_builder< FieldTCreatorT<double> >();
    }
};

}


FieldTCreator::FieldTCreator()
{
}

FieldTCreator::~FieldTCreator()
{
}


FieldTCreatorFactory::FieldTCreatorFactory(const std::string &name):
    name_(name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}


FieldTCreatorFactory::~FieldTCreatorFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}


void FieldTCreatorFactory::list(std::ostream& out) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    const char* sep = "";
    for (std::map<std::string, FieldTCreatorFactory *>::const_iterator j = m->begin() ; j != m->end() ; ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}


FieldTCreator *FieldTCreatorFactory::build(const std::string &name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, FieldTCreatorFactory *>::const_iterator j = m->find(name);

    eckit::Log::debug() << "Looking for FieldTCreatorFactory [" << name << "]" << '\n';

    if (j == m->end()) {
        eckit::Log::error() << "No FieldTCreatorFactory for [" << name << "]" << '\n';
        eckit::Log::error() << "FieldTCreatorFactories are:" << '\n';
        for (j = m->begin() ; j != m->end() ; ++j)
            eckit::Log::error() << "   " << (*j).first << '\n';
        throw eckit::SeriousBug(std::string("No FieldTCreatorFactory called ") + name);
    }

    return (*j).second->make();
}

template<> Field* FieldTCreatorT<int>::create_field( const ArrayShape& shape, const eckit::Parametrisation& params) const
  { return new FieldT<int>(shape,params); }
template<> Field* FieldTCreatorT<long>::create_field( const ArrayShape& shape, const eckit::Parametrisation& params) const
  { return new FieldT<long>(shape,params); }
template<> Field* FieldTCreatorT<float>::create_field( const ArrayShape& shape, const eckit::Parametrisation& params) const
  { return new FieldT<float>(shape,params); }
template<> Field* FieldTCreatorT<double>::create_field( const ArrayShape& shape, const eckit::Parametrisation& params) const
  { return new FieldT<double>(shape,params); }


namespace {
static FieldTCreatorBuilder< FieldTCreatorT<int> >     __FieldT__int("FieldT<int>");
static FieldTCreatorBuilder< FieldTCreatorT<long> >    __FieldT__long("FieldT<long>");
static FieldTCreatorBuilder< FieldTCreatorT<float> >   __FieldT__float("FieldT<float>");
static FieldTCreatorBuilder< FieldTCreatorT<double> >  __FieldT__double("FieldT<double>");
static FieldTCreatorBuilder< FieldTCreatorT<int> >     __FieldT__int32("FieldT<int32>");
static FieldTCreatorBuilder< FieldTCreatorT<long> >    __FieldT__int64("FieldT<int64>");
static FieldTCreatorBuilder< FieldTCreatorT<float> >   __FieldT__real32("FieldT<real32>");
static FieldTCreatorBuilder< FieldTCreatorT<double> >  __FieldT__real64("FieldT<real64>");
}

// ------------------------------------------------------------------

} // namespace fieldcreator
} // namespace atlas

