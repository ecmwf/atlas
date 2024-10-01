/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "memory_resource.h"

#include "pluto/util/Registry.h"
#include "pluto/memory_resource/MemoryPoolResource.h"

// --------------------------------------------------------------------------------------------------------

#include <cstdlib>
#include <stack>
#include "DeviceMemoryResource.h"
#include "ManagedMemoryResource.h"
#include "PinnedMemoryResource.h"
#include "pluto/host/MemoryResource.h"
#include "pluto/device/MemoryResource.h"

namespace pluto {
namespace {
struct RegisterPlutoResources {
RegisterPlutoResources() {
    register_resource("pluto::null_memory_resource", null_memory_resource());
    register_resource("pluto::new_delete_resource", new_delete_resource());
    register_resource("pluto::device_resource", device_resource());
    register_resource("pluto::managed_resource", managed_resource());
    register_resource("pluto::pinned_resource", pinned_resource());
    register_resource("pluto::pool_resource", pool_resource());
    register_resource("pluto::pinned_pool_resource", pinned_pool_resource());
    register_resource("pluto::device_pool_resource", device_pool_resource());
    register_resource("pluto::managed_pool_resource", managed_pool_resource());
}
};
}

void init() {
    static bool initialized = []() {
#if PLUTO_DEBUGGING
        std::cout << "pluto::init()" << std::endl;
#endif
        if (char* env = std::getenv("PLUTO_HOST_MEMORY_RESOURCE")) {
            pluto::host::set_default_resource(env);
        }
        if (char* env = std::getenv("PLUTO_DEVICE_MEMORY_RESOURCE")) {
            pluto::device::set_default_resource(env);
        }
        return true;
    }();
}

}

// --------------------------------------------------------------------------------------------------------

namespace pluto {

// --------------------------------------------------------------------------------------------------------

using MemoryResourceRegistry = Registry<memory_resource>;

memory_resource* register_resource(std::string_view name, memory_resource* mr) {
    return &MemoryResourceRegistry::instance().enregister(name, *mr);
}

memory_resource* register_resource(std::string_view name, std::unique_ptr<memory_resource>&& mr) {
    return &MemoryResourceRegistry::instance().enregister(name, std::move(mr));
}

void unregister_resource(std::string_view name) {
    MemoryResourceRegistry::instance().unregister(name);
}

memory_resource* get_registered_resource(std::string_view name) {
    static RegisterPlutoResources register_pluto_resources;
    return &MemoryResourceRegistry::instance().get(name);
}

bool has_registered_resource(std::string_view name) {
    return MemoryResourceRegistry::instance().has(name);
}

void unregister_resources() {
    MemoryResourceRegistry::instance().clear();
}

memory_pool_resource* pool_resource() {
    static MemoryPoolResource resource(new_delete_resource());
    return &resource;
}

// --------------------------------------------------------------------------------------------------------


struct PlutoScope {
    PlutoScope() {
        host_default_memory_resource_   = host::get_default_resource();
        device_default_memory_resource_ = device::get_default_resource();
    }
    ~PlutoScope() {
        host::set_default_resource(host_default_memory_resource_);
        device::set_default_resource(device_default_memory_resource_);
    }
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
