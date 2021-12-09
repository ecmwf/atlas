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

#include "eckit/linalg/LinearAlgebra.h"
#include "eckit/utils/Tokenizer.h"

#include "atlas/library.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

#include "atlas/runtime/Log.h"

namespace atlas {
namespace linalg {
namespace dense {

namespace {
struct backends {
    std::map<std::string, dense::Backend> map_;
    std::string current_backend_;

    static backends& instance() {
        static backends x;
        return x;
    }

    void set(const std::string& current_backend) { current_backend_ = current_backend; }

    dense::Backend& get(const std::string& type) {
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

    dense::Backend& current() { return get(current_backend_); }

private:
    backends() {
        auto configured = atlas::Library::instance().linalgDenseBackend();
        if (configured.empty()) {
            // Proposal for default backend would be "eckit_linalg":
            //     current_backend_ = configured.empty() ? backend::eckit_linalg::type() : configured;
            //
            // However to be identical in behaviour for TransLocal with atlas 0.26.0 and earlier before
            // any other codes can be adapted in short notice:
            if (eckit::linalg::LinearAlgebra::hasBackend("mkl")) {
                current_backend_ = "mkl";
            }
            else {
                current_backend_ = backend::eckit_linalg::type();
            }
        }
        map_.emplace("default", util::Config("type", current_backend_));
    }
};
}  // namespace

void current_backend(const std::string& backend) {
    backends::instance().set(backend);
}
dense::Backend& current_backend() {
    return backends::instance().current();
}

dense::Backend& default_backend(const std::string& backend) {
    return backends::instance().get(backend);
}

Backend::Backend(const std::string type): util::Config() {
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

Backend::operator std::string() const {
    return type();
}

bool Backend::available() const {
    std::string t = type();
    if (t == backend::eckit_linalg::type()) {
        if (has("backend")) {
            return eckit::linalg::LinearAlgebra::hasBackend(getString("backend"));
        }
        return true;
    }
    return eckit::linalg::LinearAlgebra::hasBackend(t);
}

}  // namespace dense
}  // namespace linalg
}  // namespace atlas
