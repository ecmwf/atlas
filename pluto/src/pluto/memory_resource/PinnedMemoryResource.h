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

#include "pluto/memory_resource.h"

namespace pluto {

// --------------------------------------------------------------------------------------------------------

class PinnedMemoryResource : public memory_resource {
public:
    PinnedMemoryResource() : PinnedMemoryResource(new_delete_resource()) {}
    PinnedMemoryResource(memory_resource* upstream) : upstream_(upstream) {}
    PinnedMemoryResource(std::string_view name) : PinnedMemoryResource(get_registered_resource(name)) {}

    void pin(void* ptr, std::size_t bytes);
    void unpin(void* ptr, std::size_t bytes);

    void* do_allocate(std::size_t bytes, std::size_t alignment) override;

    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override;
 
    bool do_is_equal(const memory_resource& other) const noexcept override;

private:
    memory_resource* upstream_;
};

memory_resource* pinned_resource();
memory_pool_resource* pinned_pool_resource();

// --------------------------------------------------------------------------------------------------------

}
