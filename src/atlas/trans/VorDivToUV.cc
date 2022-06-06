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

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/utils/Hash.h"

#include "atlas/functionspace.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/defines.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/VorDivToUV.h"

// For factory registration only
#if ATLAS_HAVE_TRANS
#include "atlas/trans/ifs/VorDivToUVIFS.h"
#define TRANS_DEFAULT "ectrans"
#else
#define TRANS_DEFAULT "local"
#endif
#include "atlas/trans/local/VorDivToUVLocal.h"  // --> recommended "local"

namespace atlas {
namespace trans {

VorDivToUVImpl::~VorDivToUVImpl() = default;

namespace {

static eckit::Mutex* local_mutex                    = nullptr;
static std::map<std::string, VorDivToUVFactory*>* m = nullptr;
static pthread_once_t once                          = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    m           = new std::map<std::string, VorDivToUVFactory*>();
}

template <typename T>
void load_builder() {
    VorDivToUVBuilder<T>("tmp");
}

struct force_link {
    force_link() {
#if ATLAS_HAVE_TRANS
        load_builder<VorDivToUVIFS>();
#endif
        load_builder<VorDivToUVLocal>();
    }
};

VorDivToUVFactory& factory(const std::string& name) {
    std::map<std::string, VorDivToUVFactory*>::const_iterator j = m->find(name);
    if (j == m->end()) {
        Log::error() << "No VorDivToUVFactory for [" << name << "]" << std::endl;
        Log::error() << "VorDivToUVFactory are:" << std::endl;
        for (j = m->begin(); j != m->end(); ++j) {
            Log::error() << "   " << (*j).first << std::endl;
        }
        throw_Exception(std::string("No VorDivToUVFactory called ") + name);
    }
    return *j->second;
}

}  // namespace

VorDivToUVFactory::VorDivToUVFactory(const std::string& name): name_(name) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ATLAS_ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}

VorDivToUVFactory::~VorDivToUVFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}

bool VorDivToUVFactory::has(const std::string& name) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    return (m->find(name) != m->end());
}

void VorDivToUVFactory::list(std::ostream& out) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    const char* sep = "";
    for (std::map<std::string, VorDivToUVFactory*>::const_iterator j = m->begin(); j != m->end(); ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}

VorDivToUV::Implementation* VorDivToUVFactory::build(const FunctionSpace& sp, const eckit::Configuration& config) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::string suffix("(" + sp.type() + ")");
    std::string name = config.getString("type", TRANS_DEFAULT) + suffix;

    Log::debug() << "Looking for TransFactory [" << name << "]" << std::endl;

    if (not config.has("type") and not has(name)) {
        name = std::string("local") + suffix;
        Log::debug() << "Looking for TransFactory [" << name << "]" << std::endl;
    }

    return factory(name).make(sp, config);
}

VorDivToUV::Implementation* VorDivToUVFactory::build(int truncation, const eckit::Configuration& config) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    static force_link static_linking;

    std::string name = config.getString("type", TRANS_DEFAULT);

    Log::debug() << "Looking for VorDivToUVFactory [" << name << "]" << std::endl;

    if (not config.has("type") and not has(name)) {
        name = std::string("local");
        Log::debug() << "Looking for VorDivToUVFactory [" << name << "]" << std::endl;
    }

    return factory(name).make(truncation, config);
}

VorDivToUV::VorDivToUV(const FunctionSpace& sp, const eckit::Configuration& config):
    Handle(VorDivToUVFactory::build(sp, config)) {}

VorDivToUV::VorDivToUV(int truncation, const eckit::Configuration& config):
    Handle(VorDivToUVFactory::build(truncation, config)) {}

int VorDivToUV::truncation() const {
    return get()->truncation();
}

// -- IFS type fields --
// These fields have special interpretation required. You need to know what
// you're doing.
// See IFS trans library.

void VorDivToUV::execute(const int nb_coeff, const int nb_fields, const double vorticity[], const double divergence[],
                         double U[], double V[], const eckit::Configuration& config) const {
    get()->execute(nb_coeff, nb_fields, vorticity, divergence, U, V, config);
}

}  // namespace trans
}  // namespace atlas
