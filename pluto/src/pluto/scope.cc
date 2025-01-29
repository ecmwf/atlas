/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "scope.h"

// --------------------------------------------------------------------------------------------------------

#include <stack>

#include "pluto/host/MemoryResource.h"
#include "pluto/device/MemoryResource.h"
#include "pluto/stream.h"

namespace pluto {

// --------------------------------------------------------------------------------------------------------

struct PlutoScope {
    PlutoScope() :
        stream_(get_current_stream()) {
        host_default_memory_resource_   = host::get_default_resource();
        device_default_memory_resource_ = device::get_default_resource();
    }
    ~PlutoScope() {
        host::set_default_resource(host_default_memory_resource_);
        device::set_default_resource(device_default_memory_resource_);
        set_current_stream(stream_);
    }
    const stream& stream_;
    memory_resource* host_default_memory_resource_;
    memory_resource* device_default_memory_resource_;
};

static std::stack<PlutoScope>& scope_stack() {
    static std::stack<PlutoScope> scope_stack_ {{PlutoScope()}};
    return scope_stack_;
}

void scope::push() {
    scope_stack().emplace();
}
void scope::pop() {
    scope_stack().pop();
}

}
