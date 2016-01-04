/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <map>
#include <string>

#include "atlas/atlas_config.h"
#include "atlas/atlas_defines.h"

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/FunctionSpace.h"
#include "atlas/numerics/Nabla.h"
#include "atlas/numerics/nabla/EdgeBasedFiniteVolume.h"

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

Nabla::Nabla(const next::FunctionSpace &fs, const eckit::Parametrisation &p)
  // : fs_(fs), config_(p)
{
}

Nabla::~Nabla()
{
}

Nabla* Nabla::create(const next::FunctionSpace &fs)
{
  return Nabla::create(fs,Config());
}

Nabla* Nabla::create(const next::FunctionSpace &fs, const eckit::Parametrisation &p)
{
  return NablaFactory::build(fs,p);
}

namespace {

template<typename T> void load_builder() { NablaBuilder<T>("tmp"); }

struct force_link {
    force_link()
    {
      load_builder< nabla::EdgeBasedFiniteVolume >();
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



Nabla* NablaFactory::build(const next::FunctionSpace& fs, const eckit::Parametrisation& p) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::map<std::string, NablaFactory *>::const_iterator j = m->find(fs.name());

    Log::debug() << "Looking for NablaFactory [" << fs.name() << "]" << '\n';

    if (j == m->end()) {
        Log::error() << "No NablaFactory for [" << fs.name() << "]" << '\n';
        Log::error() << "NablaFactories are:" << '\n';
        for (j = m->begin() ; j != m->end() ; ++j)
            Log::error() << "   " << (*j).first << '\n';
        throw eckit::SeriousBug(std::string("No NablaFactory called ") + fs.name());
    }

    return (*j).second->make(fs,p);
}

extern "C" {

void atlas__Nabla__delete(Nabla* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

Nabla* atlas__Nabla__create (const next::FunctionSpace* functionspace, const eckit::Parametrisation* params)
{
  Nabla* nabla(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(functionspace);
    ASSERT(params);
    nabla = Nabla::create(*functionspace,*params);
  );
  return nabla;
}

void atlas__Nabla__gradient (const Nabla* This, const Field* scalar, Field* grad)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(scalar);
    ASSERT(grad);
    This->gradient(*scalar,*grad);
  );
}

void atlas__Nabla__divergence (const Nabla* This, const Field* vector, Field* div)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(vector);
    ASSERT(div);
    This->divergence(*vector,*div);
  );
}

void atlas__Nabla__curl (const Nabla* This, const Field* vector, Field* curl)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(vector);
    ASSERT(curl);
    This->curl(*vector,*curl);
  );
}

void atlas__Nabla__laplacian (const Nabla* This, const Field* scalar, Field* laplacian)
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
