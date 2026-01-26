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

namespace pluto {

// --------------------------------------------------------------------------------------------------------

/// Adapts an async_memory_resource to be used by e.g. std::pmr::polymorphic_allocator
// It will always be using (de)allocate_async with the get_stream.
class AsyncMemoryResourceAdaptor : public memory_resource {
public:
    AsyncMemoryResourceAdaptor(memory_resource* mr): mr_(mr), async_mr_(dynamic_cast<async_memory_resource*>(mr)) {}

    memory_resource* upstream_resource() { return mr_; }

protected:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override;

    void do_deallocate(void* p, std::size_t bytes, std::size_t alignment) override;

    bool do_is_equal(const memory_resource& other) const noexcept override { return mr_->is_equal(other); }

private:
    memory_resource* mr_;
    async_memory_resource* async_mr_{nullptr};
};

}  // namespace pluto
