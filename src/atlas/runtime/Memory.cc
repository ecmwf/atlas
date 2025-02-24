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
#include "atlas/runtime/Exception.h"
#include "eckit/log/Bytes.h"

namespace atlas {

static bool unified_ = false;
static std::string label_;

std::string_view memory::label::get() {
    return label_;
}

void memory::label::set(std::string_view s) {
    label_ = std::string{s};
}

struct MemoryScope {
    MemoryScope() {
        pluto::scope::push();
        previous_unified_ = unified_;
    }
    ~MemoryScope() {
        unified_ = previous_unified_;
        pluto::scope::pop();
    }
    bool previous_unified_;
};

static std::stack<MemoryScope>& scope_stack() {
    static std::stack<MemoryScope> scope_stack_ {{MemoryScope()}};
    return scope_stack_;
}

void memory::set_unified(bool value) {
    unified_ = value;
}
bool memory::get_unified() {
    return unified_;
}

void memory::scope::push() {
    scope_stack().emplace();
}
void memory::scope::pop() {
    scope_stack().pop();
}


namespace memory {
context::context() {
    reset();
}

void context::reset() {
    unified_                = get_unified();
    host_memory_resource_   = pluto::host::get_default_resource();
    device_memory_resource_ = pluto::device::get_default_resource();
};

static std::map<std::string, std::unique_ptr<memory::context>> context_registry_;

bool context_exists(std::string_view name) {
    if (context_registry_.find(std::string(name)) != context_registry_.end()) {
        return true;
    }
    return false;
}

void register_context(std::string_view name) {
    std::string _name{name};
    ATLAS_ASSERT( !context_exists(name) );
    context_registry_.emplace(_name, new context());
}

void unregister_context(std::string_view name) {
    ATLAS_ASSERT( context_exists(name) );
    context_registry_.erase(std::string(name));
}

context* get_context(std::string_view name) {
    ATLAS_ASSERT( context_exists(name) );
    auto& ctx = context_registry_.at(std::string(name));
    return ctx.get();
}

void set_context(context* ctx) {
    pluto::host::set_default_resource(ctx->host_memory_resource());
    pluto::device::set_default_resource(ctx->device_memory_resource());
    set_unified(ctx->unified());
}

void set_context(std::string_view name) {
    set_context(get_context(name));
}

}

Memory::Memory(std::string_view name) :
    name_(name) {
}

void Memory::increase(size_t bytes, std::string_view label) {
    ++allocations_;
    bytes_ += bytes;
    update_maximum();
    update_largest(bytes);
    if (atlas::Library::instance().traceMemory()) {
        Log::trace() << "Memory ("<<std::setw(6) << name_<<"): " << eckit::Bytes(double(bytes_)) << "\t( +" << eckit::Bytes(double(bytes))
                        << " \t| high watermark " << eckit::Bytes(double(high_)) << "\t| label \"" << (label.empty() ? "unknown" : label ) << "\")" << std::endl;
    }
}

void Memory::decrease(size_t bytes, std::string_view label) {
    bytes_ -= bytes;
    if (atlas::Library::instance().traceMemory()) {
        Log::trace() << "Memory ("<< std::setw(6) << name_<<"): " << eckit::Bytes(double(bytes_)) << "\t( -" << eckit::Bytes(double(bytes))
                        << " \t| high watermark " << eckit::Bytes(double(high_)) << "\t| label \"" << (label.empty() ? "unknown" : label ) << "\")" << std::endl;
    }
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
 
    bool do_is_equal(const memory_resource& other) const noexcept override;

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
 
    bool do_is_equal(const memory_resource& other) const noexcept override;

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
    auto label = pluto::get_label();
    Memory::host().increase(bytes, label);
    return mr_->allocate(bytes, alignment);
}

void TraceHostMemoryResource::do_deallocate(void* p, std::size_t bytes, std::size_t alignment) {
    auto label = pluto::get_label();
    Memory::host().decrease(bytes, label);
    mr_->deallocate(p, bytes, alignment);
}

bool TraceHostMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    return mr_->is_equal(other);
}

// --------------------------------------------------------------------------------------------------------

TraceDeviceMemoryResource::TraceDeviceMemoryResource(std::string_view name, pluto::memory_resource* mr) :
    name_(name),
    mr_(mr) {
}

void* TraceDeviceMemoryResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    auto label = pluto::get_label();
    Memory::device().increase(bytes, label);
    return mr_->allocate(bytes, alignment);
}

void TraceDeviceMemoryResource::do_deallocate(void* p, std::size_t bytes, std::size_t alignment) {
    auto label = pluto::get_label();
    Memory::device().decrease(bytes, label);
    mr_->deallocate(p, bytes, alignment);
}

bool TraceDeviceMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    return mr_->is_equal(other);
}

// --------------------------------------------------------------------------------------------------------

std::unique_ptr<pluto::memory_resource> memory::host::traced_resource(pluto::memory_resource* upstream) {
    if (not upstream) {
        upstream = pluto::host::get_default_resource();
    }
    std::string name{pluto::get_registered_name(upstream)};
    return std::make_unique<pluto::TraceMemoryResource>("host   | " + name,
        std::make_unique<TraceHostMemoryResource>("host_memory   | " + name, upstream));
}

std::unique_ptr<pluto::memory_resource> memory::device::traced_resource(pluto::memory_resource* upstream) {
    if (not upstream) {
        upstream = pluto::device::get_default_resource();
    }
    std::string name{pluto::get_registered_name(upstream)};
    return std::make_unique<pluto::TraceMemoryResource>("device | " + name,
        std::make_unique<TraceDeviceMemoryResource>("device_memory | " + name, upstream));
}

// --------------------------------------------------------------------------------------------------------

}
