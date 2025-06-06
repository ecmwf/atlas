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

#include <memory>
#include <string_view>

#include "pluto/memory_resource.h"
#include "pluto/trace.h"

namespace pluto {

// --------------------------------------------------------------------------------------------------------

class TraceMemoryResource : public async_memory_resource {
public:
    TraceMemoryResource(std::string_view name, memory_resource* mr): mr_(mr), name_(name) {}

    TraceMemoryResource(memory_resource* mr): TraceMemoryResource(pluto::get_registered_name(mr), mr) {}

    TraceMemoryResource(std::string_view name): TraceMemoryResource(name, get_default_resource()) {}

    TraceMemoryResource(): TraceMemoryResource(get_default_resource()) {}

    // Take ownership of wrapped memory resource
    TraceMemoryResource(std::string_view name, std::unique_ptr<memory_resource>&& mr):
        owned_mr_(std::move(mr)), mr_(owned_mr_.get()), name_(name) {}

    // Take ownership of wrapped memory resource
    TraceMemoryResource(std::unique_ptr<memory_resource>&& mr):
        owned_mr_(std::move(mr)), mr_(owned_mr_.get()), name_(get_registered_name(mr_)) {}

    memory_resource* upstream_resource() { return mr_; }

protected:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override;

    void do_deallocate(void* p, std::size_t bytes, std::size_t alignment) override;

    void* do_allocate_async(std::size_t bytes, std::size_t alignment, stream_view) override;

    void do_deallocate_async(void* p, std::size_t bytes, std::size_t alignment, stream_view) override;

    bool do_is_equal(const memory_resource& other) const noexcept override { return mr_->is_equal(other); }

private:
    std::unique_ptr<memory_resource> owned_mr_;
    memory_resource* mr_;
    std::string name_;
    static int nest;
};

}  // namespace pluto
