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

#include "pluto/memory_resource/memory_resource.h"

//------------------------------------------------------------------------------

namespace atlas {

// --------------------------------------------------------------------------------------------------------

class Memory {
public:
    static Memory& instance() {
        return host();
    }

    static Memory& host() {
        static Memory _instance("host");
        return _instance;
    }

    static Memory& device() {
        static Memory _instance("device");
        return _instance;
    }

    Memory& operator+=(size_t bytes);

    Memory& operator-=(size_t bytes);

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

pluto::memory_resource* host_memory_resource(pluto::memory_resource* upstream = nullptr);
pluto::memory_resource* device_memory_resource(pluto::memory_resource* upstream = nullptr);

// --------------------------------------------------------------------------------------------------------

}  // namespace atlas
