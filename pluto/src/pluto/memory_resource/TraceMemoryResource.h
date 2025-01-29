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

#include <string_view>
#include <memory>

#include "pluto/memory_resource.h"
#include "pluto/trace.h"

namespace pluto {

// --------------------------------------------------------------------------------------------------------

class TraceMemoryResource : public async_memory_resource {
public:

    TraceMemoryResource( std::string_view name, memory_resource* mr ) :
        name_(name),
        mr_(mr) {
    }

    TraceMemoryResource( std::string_view name) :
        name_(name),
        mr_(get_default_resource()) {
    }

    // Take ownership of wrapped memory resource
    TraceMemoryResource( std::string_view name, std::unique_ptr<memory_resource>&& mr ) :
        name_(name),
        owned_mr_(std::move(mr)),
        mr_(owned_mr_.get()) {
    }

    memory_resource* upstream_resource() {
        return mr_;
    }

protected:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override;
 
    void do_deallocate(void* p, std::size_t bytes, std::size_t alignment) override;
 
    void* do_allocate_async(std::size_t bytes, std::size_t alignment, const stream&) override;
 
    void do_deallocate_async(void* p, std::size_t bytes, std::size_t alignment, const stream&) override;

    bool do_is_equal(const memory_resource& other) const noexcept override {
        return mr_->is_equal(other);
    }

private:
    std::string name_;
    std::unique_ptr<memory_resource> owned_mr_;
    memory_resource* mr_;
    static int nest;
};

}
