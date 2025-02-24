/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <atomic>
#include <string>
#include <string_view>

#include "pluto/pluto.h"

//------------------------------------------------------------------------------

namespace atlas {

// --------------------------------------------------------------------------------------------------------

class Memory {
public:
    static Memory& instance() {
        return host();
    }

    static Memory& host() {
        static Memory _instance(" host ");
        return _instance;
    }

    static Memory& device() {
        static Memory _instance("device");
        return _instance;
    }

    void increase(size_t bytes, std::string_view label);
    void decrease(size_t bytes, std::string_view label);

    // Memory& operator+=(size_t bytes);

    // Memory& operator-=(size_t bytes);

    size_t allocated() const {
        return bytes_;
    }

    size_t highWatermark() const {
        return high_;
    }

    /// Number of allocations
    size_t allocations() const {
        return allocations_;
    }

    size_t largestAllocation() const {
        return largest_;
    }

    void reset();

private:
    Memory(std::string_view name);
    void update_maximum() noexcept;
    void update_largest(size_t bytes) noexcept;

    std::atomic<size_t> allocations_{0};
    std::atomic<size_t> bytes_{0};
    std::atomic<size_t> high_{0};
    std::atomic<size_t> largest_{0};
    std::string name_;
};

// --------------------------------------------------------------------------------------------------------

namespace memory {

namespace host {
// inline void set_default_resource(pluto::memory_resource* mr) {
//     pluto::host::set_default_resource(mr);
// }

// inline void set_default_resource(std::string_view name) {
//     pluto::host::set_default_resource(name);
// }

std::unique_ptr<pluto::memory_resource> traced_resource(pluto::memory_resource* upstream = nullptr);
}

namespace device {
// inline void set_default_resource(pluto::memory_resource* mr) {
//     pluto::device::set_default_resource(mr);
// }

// inline void set_default_resource(std::string_view name) {
//     pluto::device::set_default_resource(name);
// }

std::unique_ptr<pluto::memory_resource> traced_resource(pluto::memory_resource* upstream = nullptr);
}

bool get_unified();
void set_unified(bool);

struct scope {
    scope() {
        push();
    }
    ~scope() {
        pop();
    }
    static void push();
    static void pop();
};

class context {
public:
    context();
    pluto::memory_resource* host_memory_resource() { return host_memory_resource_; }
    pluto::memory_resource* device_memory_resource() { return device_memory_resource_; }
    bool unified() { return unified_; }
    void reset();
private:
    bool unified_;
    pluto::memory_resource* host_memory_resource_;
    pluto::memory_resource* device_memory_resource_;
};

void register_context(std::string_view name);
void unregister_context(std::string_view name);
bool context_exists(std::string_view name);
void set_context(std::string_view name);
void set_context(context* ctx);
context* get_context(std::string_view name);

class label {
public:
    label(std::string_view s) {
        previous_ = get();
        set(s);
    }
    ~label() {
        set(previous_);
    }
    static std::string_view get();
    static void set(std::string_view);
private:
    std::string previous_;
};


}

}  // namespace atlas
