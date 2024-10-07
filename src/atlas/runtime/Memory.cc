/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Memory.h"

#include <iostream>
#include <memory>
#include <iomanip>
#include <stack>

#include "pluto/pluto.h"

#include "atlas/library/Library.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Log.h"
#include "eckit/log/Bytes.h"

namespace atlas {

struct MemoryScope {
    MemoryScope() {
        pluto::scope::push();
    }
    ~MemoryScope() {
        pluto::scope::pop();
    }
    MemoryScope(const MemoryScope& previous) {
        device_memory_mapped_ = previous.device_memory_mapped_;
    }
    bool device_memory_mapped_ = false;
};

static std::stack<MemoryScope>& scope_stack() {
    static std::stack<MemoryScope> scope_stack_ {{MemoryScope()}};
    return scope_stack_;
}

void memory::set_unified(bool value) {
    scope_stack().top().device_memory_mapped_ = value;
}
bool memory::get_unified() {
    return scope_stack().top().device_memory_mapped_;
}

void memory::scope::push() {
    scope_stack().emplace(scope_stack().top());
}
void memory::scope::pop() {
    scope_stack().pop();
}

Memory::Memory(std::string_view name) :
    name_(name) {
}

Memory& Memory::operator+=(size_t bytes) {
    ++allocations_;
    bytes_ += bytes;
    update_maximum();
    update_largest(bytes);
    if (atlas::Library::instance().traceMemory()) {
        Log::trace() << "Memory ("<<std::setw(6) << name_<<"): " << eckit::Bytes(double(bytes_)) << "\t( +" << eckit::Bytes(double(bytes))
                        << " \t| high watermark " << eckit::Bytes(double(high_)) << "\t)" << std::endl;
    }
    return *this;
}

Memory& Memory::operator-=(size_t bytes) {
    bytes_ -= bytes;
    if (atlas::Library::instance().traceMemory()) {
        Log::trace() << "Memory ("<< std::setw(6) << name_<<"): " << eckit::Bytes(double(bytes_)) << "\t( -" << eckit::Bytes(double(bytes))
                        << " \t| high watermark " << eckit::Bytes(double(high_)) << "\t)" << std::endl;
    }
    return *this;
}

void Memory::update_maximum() noexcept {
    size_t prev_value = high_;
    while (prev_value < bytes_ && !high_.compare_exchange_weak(prev_value, bytes_)) {
    }
}

void Memory::update_largest(std::size_t largest) noexcept {
    size_t prev_value = largest_;
    while (prev_value < largest && !largest_.compare_exchange_weak(prev_value, largest)) {
    }
}


// --------------------------------------------------------------------------------------------------------

class TraceHostMemoryResource : public pluto::memory_resource {
public:
    TraceHostMemoryResource(std::string_view name, pluto::memory_resource* mr);

    pluto::memory_resource* upstream() { return mr_; }
protected:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override;

    void do_deallocate(void* p, std::size_t bytes, std::size_t alignment) override;
 
    bool do_is_equal(const pluto::memory_resource& other) const noexcept override;

private:
    std::string name_;
    pluto::memory_resource* mr_;
};

// --------------------------------------------------------------------------------------------------------

class TraceDeviceMemoryResource : public pluto::memory_resource {
public:
    TraceDeviceMemoryResource(std::string_view name, pluto::memory_resource* mr);

    pluto::memory_resource* upstream() { return mr_; }
protected:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override;

    void do_deallocate(void* p, std::size_t bytes, std::size_t alignment) override;
 
    bool do_is_equal(const pluto::memory_resource& other) const noexcept override;

private:
    std::string name_;
    pluto::memory_resource* mr_;
};

// --------------------------------------------------------------------------------------------------------

TraceHostMemoryResource::TraceHostMemoryResource(std::string_view name, pluto::memory_resource* mr) :
    name_(name),
    mr_(mr) {
}

void* TraceHostMemoryResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    Memory::host() += bytes;
    return mr_->allocate(bytes, alignment);
}

void TraceHostMemoryResource::do_deallocate(void* p, std::size_t bytes, std::size_t alignment) {
    Memory::host() -= bytes;
    mr_->deallocate(p, bytes, alignment);
}

bool TraceHostMemoryResource::do_is_equal(const pluto::memory_resource& other) const noexcept {
    return mr_->is_equal(other);
}

// --------------------------------------------------------------------------------------------------------

TraceDeviceMemoryResource::TraceDeviceMemoryResource(std::string_view name, pluto::memory_resource* mr) :
    name_(name),
    mr_(mr) {
}

void* TraceDeviceMemoryResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    Memory::device() += bytes;
    return mr_->allocate(bytes, alignment);
}

void TraceDeviceMemoryResource::do_deallocate(void* p, std::size_t bytes, std::size_t alignment) {
    Memory::device() -= bytes;
    mr_->deallocate(p, bytes, alignment);
}

bool TraceDeviceMemoryResource::do_is_equal(const pluto::memory_resource& other) const noexcept {
    return mr_->is_equal(other);
}

// --------------------------------------------------------------------------------------------------------

std::unique_ptr<pluto::memory_resource> memory::host::traced_resource(pluto::memory_resource* upstream) {
    if (not upstream) {
        upstream = pluto::host::get_default_resource();
    } 
    return std::make_unique<pluto::TraceMemoryResource>("host_memory", 
        std::make_unique<TraceHostMemoryResource>("host_memory", upstream));
}

std::unique_ptr<pluto::memory_resource> memory::device::traced_resource(pluto::memory_resource* upstream) {
    if (not upstream) {
        upstream = pluto::device::get_default_resource();
    } 
    return std::make_unique<pluto::TraceMemoryResource>("device_memory", 
        std::make_unique<TraceDeviceMemoryResource>("device_memory", upstream));
}

// --------------------------------------------------------------------------------------------------------

}
