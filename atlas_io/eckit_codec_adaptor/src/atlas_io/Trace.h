/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once


#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <limits>

namespace eckit {
class CodeLocation;
}

namespace atlas {
namespace io {

struct TraceHook {
    TraceHook()          = default;
    virtual ~TraceHook() = default;
};

struct TraceHookRegistry {
    using TraceHookBuilder = std::function<std::unique_ptr<TraceHook>(const eckit::CodeLocation&, const std::string&)>;
    std::vector<TraceHookBuilder> hooks;
    std::vector<int> enabled_;
    static TraceHookRegistry& instance() {
        static TraceHookRegistry instance;
        return instance;
    }
    static size_t add(TraceHookBuilder&& hook) {
        instance().hooks.emplace_back(hook);
        instance().enabled_.emplace_back(true);
        return instance().hooks.size() - 1;
    }
    static size_t add(const TraceHookBuilder& hook) {
        instance().hooks.emplace_back(hook);
        instance().enabled_.emplace_back(true);
        return instance().hooks.size() - 1;
    }
    static void enable(size_t id) { instance().enabled_[id] = true; }
    static void disable(size_t id) { instance().enabled_[id] = false; }
    static bool enabled(size_t id) { return instance().enabled_[id]; }
    static size_t size() { return instance().hooks.size(); }
    static TraceHookBuilder& hook(size_t id) { return instance().hooks[id]; }
    static size_t invalidId() { return std::numeric_limits<size_t>::max(); }

private:
    TraceHookRegistry() = default;
};

struct Trace {
    using Labels = std::vector<std::string>;
    Trace(const eckit::CodeLocation& loc);
    Trace(const eckit::CodeLocation& loc, const std::string& title);
    Trace(const eckit::CodeLocation& loc, const std::string& title, const Labels& labels);

private:
    std::vector<std::unique_ptr<TraceHook>> hooks_;
};

}  // namespace io
}  // namespace atlas

#include "atlas_io/detail/BlackMagic.h"

#define ATLAS_IO_TRACE(...) __ATLAS_IO_TYPE(::atlas::io::Trace, Here() __ATLAS_IO_COMMA_ARGS(__VA_ARGS__))
#define ATLAS_IO_TRACE_SCOPE(...) __ATLAS_IO_TYPE_SCOPE(::atlas::io::Trace, Here() __ATLAS_IO_COMMA_ARGS(__VA_ARGS__))
