/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <map>
#include <mutex>
#include <ostream>
#include <string>
#include <vector>

#include "TransFactory.h"

#include "atlas/functionspace.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/Cache.h"
#include "atlas/trans/Trans.h"

#include "atlas/trans/local/TransLocal.h"
#if ATLAS_HAVE_TRANS
#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/ifs/TransIFSNodeColumns.h"
#include "atlas/trans/ifs/TransIFSStructuredColumns.h"
#endif

namespace atlas {
namespace trans {

//----------------------------------------------------------------------------------------------------------------------

void force_link() {
    static struct Link {
        Link() {
            TransBuilderGrid<TransLocal>();
#if ATLAS_HAVE_TRANS
            TransBuilderGrid<TransIFS>();
            TransBuilderFunctionSpace<TransIFSStructuredColumns>();
            TransBuilderFunctionSpace<TransIFSNodeColumns>();
#endif
        }
    } link;
}

//----------------------------------------------------------------------------------------------------------------------

namespace {
struct default_backend {
#if ATLAS_HAVE_TRANS
    std::string value = "ectrans";
#else
    std::string value = "local";
#endif
    static default_backend instance() {
        static default_backend x;
        return x;
    }

private:
    default_backend() = default;
};
}  // namespace

class TransBackend {
    using lock_guard = std::lock_guard<std::mutex>;

protected:
    TransBackend() { default_options_ = util::Config("type", default_backend::instance().value); }

private:
    mutable std::mutex mutex_;
    std::map<std::string, int> backends_;

    util::Config default_options_;

public:
    std::vector<std::string> keys() const {
        lock_guard lock(mutex_);
        std::vector<std::string> _keys;
        _keys.reserve(backends_.size());
        for (const auto& key_value : backends_) {
            _keys.emplace_back(key_value.first);
        }
        return _keys;
    }
    void list(std::ostream& out) const {
        lock_guard lock(mutex_);
        const char* sep = "";
        for (const auto& map_pair : backends_) {
            out << sep << map_pair.first;
            sep = ", ";
        }
    }
    bool has(const std::string& backend) const {
        lock_guard lock(mutex_);
        return (backends_.find(backend) != backends_.end());
    }
    void add(const std::string& backend) {
        lock_guard lock(mutex_);
        if (backends_.find(backend) == backends_.end()) {
            backends_[backend] = 1;
        }
        else {
            ++backends_[backend];
        }
    }
    void remove(const std::string& backend) {
        lock_guard lock(mutex_);
        --backends_[backend];
        if (backends_[backend] == 0) {
            backends_.erase(backend);
        }
    }

    void backend(const std::string& backend) {
        lock_guard lock(mutex_);
        default_options_.set("type", backend);
    }
    std::string backend() { return default_options_.getString("type"); }
    void config(const eckit::Configuration& config) {
        std::string type = default_options_.getString("type");
        default_options_ = config;
        if (not config.has("type")) {
            default_options_.set("type", type);
        }
    }
    const eckit::Configuration& config() { return default_options_; }


public:
    static TransBackend& instance() {
        static TransBackend env;
        return env;
    }
};


TransFactory::TransFactory(const std::string& name, const std::string& backend):
    Factory(name), name_(name), backend_(backend) {
    TransBackend::instance().add(backend);
}

TransFactory::~TransFactory() {
    TransBackend::instance().remove(backend_);
}

void TransFactory::list(std::ostream& out) {
    TransBackend::instance().list(out);
}

bool TransFactory::has(const std::string& backend) {
    return TransBackend::instance().has(backend);
}

void TransFactory::backend(const std::string& backend) {
    TransBackend::instance().backend(backend);
}

std::string TransFactory::backend() {
    return TransBackend::instance().backend();
}

const eckit::Configuration& TransFactory::config() {
    return TransBackend::instance().config();
}

void TransFactory::config(const eckit::Configuration& config) {
    TransBackend::instance().config(config);
}

const TransImpl* TransFactory::build(const FunctionSpace& gp, const FunctionSpace& sp,
                                     const eckit::Configuration& config) {
    return build(Cache(), gp, sp, config);
}

const TransImpl* TransFactory::build(const Cache& cache, const FunctionSpace& gp, const FunctionSpace& sp,
                                     const eckit::Configuration& config) {
    force_link();

    if (cache.trans()) {
        Log::debug() << "Creating Trans from cache, ignoring any other arguments" << std::endl;
        return cache.trans();
    }

    util::Config options = TransBackend::instance().config();
    options.set(eckit::LocalConfiguration(config));

    std::string backend = options.getString("type");

    Log::debug() << "Looking for TransFactory [" << backend << "]" << std::endl;
    if (!TransBackend::instance().has(backend)) {
        Log::error() << "No TransFactory for [" << backend << "]" << std::endl;
        Log::error() << "TransFactories are :" << std::endl;
        TransBackend::instance().list(Log::error());
        Log::error() << std::endl;
        throw_Exception(std::string("No TransFactory called ") + backend);
    }

    std::string suffix("(" + gp.type() + "," + sp.type() + ")");
    std::string builder = backend + suffix;

    Log::debug() << "Looking for TransFactory [" << builder << "]" << std::endl;
    auto factory = get(builder);
    return factory->make(cache, gp, sp, options);
}

const TransImpl* TransFactory::build(const Grid& grid, int truncation, const eckit::Configuration& config) {
    return build(Cache(), grid, truncation, config);
}

const TransImpl* TransFactory::build(const Grid& grid, const Domain& domain, int truncation,
                                     const eckit::Configuration& config) {
    return build(Cache(), grid, domain, truncation, config);
}

const TransImpl* TransFactory::build(const Cache& cache, const Grid& grid, int truncation,
                                     const eckit::Configuration& config) {
    return build(cache, grid, grid.domain(), truncation, config);
}

const TransImpl* TransFactory::build(const Cache& cache, const Grid& grid, const Domain& domain, int truncation,
                                     const eckit::Configuration& config) {
    force_link();

    if (cache.trans()) {
        Log::debug() << "Creating Trans from cache, ignoring any other arguments" << std::endl;
        return cache.trans();
    }
    util::Config options = TransBackend::instance().config();
    options.set(eckit::LocalConfiguration(config));

    std::string backend = options.getString("type");

    Log::debug() << "Looking for TransFactory [" << backend << "]" << std::endl;
    if (!TransBackend::instance().has(backend)) {
        Log::error() << "No TransFactory for [" << backend << "]" << std::endl;
        Log::error() << "TransFactories are :" << std::endl;
        TransBackend::instance().list(Log::error());
        Log::error() << std::endl;
        throw_Exception(std::string("No TransFactory called ") + backend);
    }
    std::string builder = backend;
    auto factory        = get(builder);
    return factory->make(cache, grid, domain, truncation, options);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
