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
#include "atlas/internals/atlas_config.h"
#include "atlas/internals/atlas_defines.h"
#include "atlas/numerics/Nabla.h"
#include "atlas/numerics/Method.h"
#include "atlas/numerics/fvm/Nabla.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"

namespace {

    static eckit::Mutex *local_mutex = 0;
    static std::map<std::string, atlas::numerics::NablaFactory *> *m = 0;
    static pthread_once_t once = PTHREAD_ONCE_INIT;

    static void init() {
        local_mutex = new eckit::Mutex();
        m = new std::map<std::string, atlas::numerics::NablaFactory *>();
    }
}


namespace atlas {
namespace numerics {

Nabla::Nabla(const Method &method, const eckit::Parametrisation &p)
{
}

Nabla::~Nabla()
{
}

Nabla* Nabla::create(const Method &method)
{
  return Nabla::create(method,util::Config());
}

Nabla* Nabla::create(const Method &method, const eckit::Parametrisation &p)
{
  return NablaFactory::build(method,p);
}

namespace {

template<typename T> void load_builder() { NablaBuilder<T>("tmp"); }

struct force_link {
    force_link()
    {
      load_builder< fvm::Nabla >();
    }
};

}


NablaFactory::NablaFactory(const std::string &name):
    name_(name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;

}


NablaFactory::~NablaFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}


void NablaFactory::list(std::ostream& out) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    const char* sep = "";
    for (std::map<std::string, NablaFactory *>::const_iterator j = m->begin() ; j != m->end() ; ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}

bool NablaFactory::has(const std::string& name)
{
  pthread_once(&once, init);

  eckit::AutoLock<eckit::Mutex> lock(local_mutex);

  static force_link static_linking;

  return ( m->find(name) != m->end() );
}



Nabla* NablaFactory::build(const Method& method, const eckit::Parametrisation& p) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, NablaFactory *>::const_iterator j = m->find(method.name());

    Log::debug() << "Looking for NablaFactory [" << method.name() << "]" << '\n';

    if (j == m->end()) {
        Log::error() << "No NablaFactory for [" << method.name() << "]" << '\n';
        Log::error() << "NablaFactories are:" << '\n';
        for (j = m->begin() ; j != m->end() ; ++j)
            Log::error() << "   " << (*j).first << '\n';
        throw eckit::SeriousBug(std::string("No NablaFactory called ") + method.name());
    }

    return (*j).second->make(method,p);
}

extern "C" {

void atlas__Nabla__delete(Nabla* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

Nabla* atlas__Nabla__create (const Method* method, const eckit::Parametrisation* params)
{
  Nabla* nabla(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(method);
    ASSERT(params);
    nabla = Nabla::create(*method,*params);
  );
  return nabla;
}

void atlas__Nabla__gradient (const Nabla* This, const field::Field* scalar, field::Field* grad)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(scalar);
    ASSERT(grad);
    This->gradient(*scalar,*grad);
  );
}

void atlas__Nabla__divergence (const Nabla* This, const field::Field* vector, field::Field* div)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(vector);
    ASSERT(div);
    This->divergence(*vector,*div);
  );
}

void atlas__Nabla__curl (const Nabla* This, const field::Field* vector, field::Field* curl)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(vector);
    ASSERT(curl);
    This->curl(*vector,*curl);
  );
}

void atlas__Nabla__laplacian (const Nabla* This, const field::Field* scalar, field::Field* laplacian)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(scalar);
    ASSERT(laplacian);
    This->laplacian(*scalar,*laplacian);
  );
}


}


} // namespace numerics
} // namespace atlas
