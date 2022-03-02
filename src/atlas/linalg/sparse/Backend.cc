/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Backend.h"

#include <map>

#include "atlas/library/config.h"
#if ATLAS_ECKIT_HAVE_ECKIT_585
#include "eckit/linalg/LinearAlgebraSparse.h"
#else
#include "eckit/linalg/LinearAlgebra.h"
#endif
#include "eckit/utils/Tokenizer.h"

#include "atlas/library.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace linalg {
namespace sparse {

namespace {


struct backends {
    std::map<std::string, sparse::Backend> map_;
    std::string current_backend_;

    static backends& instance() {
        static backends x;
        return x;
    }

    void set(const std::string& current_backend) { current_backend_ = current_backend; }

    sparse::Backend& get(const std::string& type) {
        if (map_.find(type) == map_.end()) {
            std::vector<std::string> tokens;
            eckit::Tokenizer{'.'}(type, tokens);
            ATLAS_ASSERT(tokens.size() <= 2);
            if (tokens.size() == 1) {
                map_.emplace(type, util::Config("type", type));
            }
            else {
                util::Config b;
                b.set("type", tokens[0]);
                b.set("backend", tokens[1]);
                map_.emplace(type, b);
            }
        }
        return map_[type];
    }

    sparse::Backend& get() { return get(current_backend_); }

private:
    backends() {
        auto configured  = atlas::Library::instance().linalgSparseBackend();
        current_backend_ = configured.empty() ? backend::openmp::type() : configured;
    }
};
}  // namespace

void current_backend(const std::string& backend) {
    backends::instance().set(backend);
}
sparse::Backend& current_backend() {
    return backends::instance().get();
}

sparse::Backend& default_backend(const std::string& backend) {
    return backends::instance().get(backend);
}

Backend::Backend(const std::string& type): util::Config() {
    if (not type.empty()) {
        set(default_backend(type));
    }
    else {
        set(current_backend());
    }
}

Backend::Backend(const eckit::Configuration& other): util::Config(other) {
    ATLAS_ASSERT(has("type"));
}

bool Backend::available() const {
    std::string t = type();
    if (t == backend::openmp::type()) {
        return true;
    }
    if (t == backend::eckit_linalg::type()) {
        if (has("backend")) {
#if ATLAS_ECKIT_HAVE_ECKIT_585
            return eckit::linalg::LinearAlgebraSparse::hasBackend(getString("backend"));
#else
            return eckit::linalg::LinearAlgebra::hasBackend(getString("backend"));
#endif
        }
        return true;
    }
#if ATLAS_ECKIT_HAVE_ECKIT_585
    return eckit::linalg::LinearAlgebraSparse::hasBackend(t);
#else
    return eckit::linalg::LinearAlgebra::hasBackend(t);
#endif
}


}  // namespace sparse
}  // namespace linalg
}  // namespace atlas
