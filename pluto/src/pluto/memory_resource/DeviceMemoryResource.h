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

#include "pluto/memory_resource/memory_resource.h"
#include "pluto/util/Alignment.h"

namespace pluto {

// --------------------------------------------------------------------------------------------------------

class DeviceMemoryResource : public memory_resource {
public:
    using alignment_t = std::size_t;
    static constexpr alignment_t alignment = default_alignment();

    DeviceMemoryResource(){};

    ~DeviceMemoryResource(){};

public:

    void* do_allocate(std::size_t bytes, alignment_t) override;
    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override;
    bool do_is_equal(const memory_resource& other) const noexcept override;

};

memory_resource* device_resource();
memory_pool_resource* device_pool_resource();

// --------------------------------------------------------------------------------------------------------

}
