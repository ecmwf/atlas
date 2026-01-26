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

#include <cstddef>
#include <string>
#include <string_view>

namespace pluto {
class memory_tracker {
public:
    memory_tracker(std::string_view name);
    ~memory_tracker();

    void allocate(std::size_t bytes);

    void deallocate(std::size_t bytes);

    std::size_t total_allocations() const;

    std::size_t allocations() const;

    std::size_t bytes() const;

    // High watermark of allocation in bytes
    std::size_t high_watermark() const;

    // Largest individual allocation in bytes
    std::size_t largest_allocation() const;

    void reset();

    bool enable(bool value) {
        bool previously_enabled = enabled_;
        enabled_ = value;
        return previously_enabled;
    }

    std::string_view name() const;

private:
    std::string name_;
    bool enabled_{true};
    struct memory_tracker_data;
    memory_tracker_data* data_{nullptr};
};

} // namespace pluto

namespace pluto::memory {

extern memory_tracker host;
extern memory_tracker host_pool;
extern memory_tracker pinned;
extern memory_tracker pinned_pool;
extern memory_tracker device;
extern memory_tracker device_pool;
extern memory_tracker managed;
extern memory_tracker managed_pool;

std::string report();
std::string report(std::string_view prefix);

} // namespace pluto::memory

