/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "memory.h"

#include <atomic>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <vector>

#include "pluto.h"

namespace pluto {

struct memory_tracker::memory_tracker_data {
    void update_high() noexcept{
        std::size_t prev_value = high_;
        while (prev_value < bytes_ && !high_.compare_exchange_weak(prev_value, bytes_)) {}
    }
    void update_largest(std::size_t largest) noexcept {
        size_t prev_value = largest_;
        while (prev_value < largest && !largest_.compare_exchange_weak(prev_value, largest)) {}
    }
    std::atomic<std::size_t> total_allocations_{0};
    std::atomic<std::size_t> allocations_{0};
    std::atomic<std::size_t> bytes_{0};
    std::atomic<std::size_t> high_{0};
    std::atomic<std::size_t> largest_{0};
};

memory_tracker::memory_tracker(std::string_view name) : name_(name) {
    data_ = new memory_tracker::memory_tracker_data();
}

memory_tracker::~memory_tracker() {
    delete data_;
}

void memory_tracker::allocate(std::size_t bytes) {
    if (enabled_) {
        ++(data_->allocations_);
        ++(data_->total_allocations_);
        data_->bytes_ += bytes;
        data_->update_high();
        data_->update_largest(bytes);
    }
}

void memory_tracker::deallocate(std::size_t bytes) {
    if (enabled_) {
        data_->bytes_ -= bytes;
        --(data_->allocations_);
    }
}

std::size_t memory_tracker::bytes() const {
    return data_->bytes_;
}

std::string_view memory_tracker::name() const {
    return name_;
}


std::size_t memory_tracker::high_watermark() const {
    return data_->high_;
}

std::size_t memory_tracker::total_allocations() const {
    return data_->total_allocations_;
}

std::size_t memory_tracker::allocations() const {
    return data_->allocations_;
}

std::size_t memory_tracker::largest_allocation() const {
    return data_->largest_;
}

void memory_tracker::reset() {
    data_->total_allocations_ = 0;
    data_->allocations_ = 0;
    data_->bytes_ = 0;
    data_->high_ = 0;
    data_->largest_ = 0;
}

}

namespace pluto::memory {

memory_tracker host("host");
memory_tracker host_pool("host_pool");
memory_tracker pinned("pinned");
memory_tracker pinned_pool("pinned_pool");
memory_tracker device("device");
memory_tracker device_pool("device_pool");
memory_tracker managed("managed");
memory_tracker managed_pool("managed_pool");


std::string report(std::string_view prefix) {

    using trace::format_bytes;

    auto format_allocation = [](std::size_t size, std::size_t bytes) -> std::string {
        std::stringstream ss;
        ss << format_bytes(bytes) << " (" << size << ")";
        std::string str = ss.str();
        return str;
    };

    std::vector<memory_tracker*> memory_trackers {
        &host,
        &host_pool,
        &pinned,
        &pinned_pool,
        &device,
        &device_pool,
        &managed,
        &managed_pool,
    };

    int first_colwidth = 25;
    int colwidth = 15;

    std::stringstream out;
    out << prefix << std::setw(first_colwidth) << std::left << "Memory resource";
    for( auto* mem : memory_trackers) {
        out << std::setw(colwidth) << mem->name();
    }
    out << '\n';

    out << prefix << std::setw(first_colwidth) << std::left << "Total allocations:";
    for( auto* mem : memory_trackers) {
        out << std::setw(colwidth) << mem->total_allocations();
    }
    out << '\n';


    out << prefix << std::setw(first_colwidth) << std::left << "High Watermark:";
    for( auto* mem : memory_trackers) {
        out << std::setw(colwidth) << format_bytes(mem->high_watermark());
    }
    out << '\n';

    out << prefix << std::setw(first_colwidth) << std::left << "Largest Allocation:";
    for( auto* mem : memory_trackers) {
        out << std::setw(colwidth) << format_bytes(mem->largest_allocation());
    }
    out << '\n';

    out << prefix << std::setw(first_colwidth) << std::left << "Allocated:";
    for( auto* mem : memory_trackers) {
        out << std::setw(colwidth) << format_allocation(mem->allocations(), mem->bytes());
    }
    out << '\n';

    out << prefix << std::setw(first_colwidth) << std::left << "Capacity:";
    for( auto* mem : memory_trackers) {
        if (mem->name().find("host_pool") != std::string::npos) {
            out << std::setw(colwidth) << format_bytes(pluto::host_pool_resource()->capacity());
        }
        else if (mem->name().find("pinned_pool") != std::string::npos) {
            out << std::setw(colwidth) << format_bytes(pluto::pinned_pool_resource()->capacity());
        }
        else if (mem->name().find("device_pool") != std::string::npos) {
            out << std::setw(colwidth) << format_bytes(pluto::device_pool_resource()->capacity());
        }
        else if (mem->name().find("managed_pool") != std::string::npos) {
            out << std::setw(colwidth) << format_bytes(pluto::managed_pool_resource()->capacity());
        }
        else {
            out << std::setw(colwidth) << "-";
        }
    }
    out << '\n';

    return out.str();
}

std::string report() {
    return report("");
}

}  // namespace pluto::memory

