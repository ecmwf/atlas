/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Trace.h"

#include "eckit/log/CodeLocation.h"

namespace atlas {
namespace io {

atlas::io::Trace::Trace(const eckit::CodeLocation& loc) {
    for (auto& hook : TraceHookRegistry::instance().hooks) {
        hooks_.emplace_back(hook(loc, loc.func()));
    }
}

Trace::Trace(const eckit::CodeLocation& loc, const std::string& title) {
    for (auto& hook : TraceHookRegistry::instance().hooks) {
        hooks_.emplace_back(hook(loc, title));
    }
}

Trace::Trace(const eckit::CodeLocation& loc, const std::string& title, const Labels& labels) {
    for (auto& hook : TraceHookRegistry::instance().hooks) {
        hooks_.emplace_back(hook(loc, title));
    }
}

}  // namespace io
}  // namespace atlas
