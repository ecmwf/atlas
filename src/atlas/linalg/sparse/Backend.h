/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>
#include "atlas/util/Config.h"

namespace atlas {
namespace linalg {
namespace sparse {

struct Backend;

void current_backend(const std::string& backend);
sparse::Backend& current_backend();
sparse::Backend& default_backend(const std::string& backend);


struct Backend : util::Config {
    Backend(): util::Config() { set(current_backend()); }
    Backend(const std::string& type);
    Backend(const eckit::Configuration& other);
    std::string type() const { return getString("type"); }
    operator std::string() const { return type(); }
    bool available() const;
};

namespace backend {
struct openmp : Backend {
    static std::string type() { return "openmp"; }
    openmp(): Backend(type()) {}
};

struct eckit_linalg : Backend {
    static std::string type() { return "eckit_linalg"; }
    eckit_linalg(): Backend(type()) {}
};

struct hicsparse : Backend {
    static std::string type() { return "hicsparse"; }
    hicsparse(): Backend(type()) {}
};
}  // namespace backend


}  // namespace sparse
}  // namespace linalg
}  // namespace atlas
