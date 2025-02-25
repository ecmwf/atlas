/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "pluto/memory_resource.h"

#include "pluto/detail/Registry.h"
#include "pluto/memory_resource/MemoryPoolResource.h"

// --------------------------------------------------------------------------------------------------------

#include <cstdlib>
#include "DeviceMemoryResource.h"
#include "ManagedMemoryResource.h"
#include "PinnedMemoryResource.h"
#include "pluto/device/MemoryResource.h"
#include "pluto/host/MemoryResource.h"

namespace pluto {
namespace {
struct RegisterPlutoResources {
private:
    RegisterPlutoResources() {
        register_resource("pluto::null_memory_resource", null_memory_resource());
        register_resource("pluto::new_delete_resource", new_delete_resource());
        register_resource("pluto::pinned_resource", pinned_resource());
        register_resource("pluto::device_resource", device_resource());
        register_resource("pluto::managed_resource", managed_resource());
        register_resource("pluto::pinned_pool_resource", pinned_pool_resource());
        register_resource("pluto::device_pool_resource", device_pool_resource());
        register_resource("pluto::managed_pool_resource", managed_pool_resource());
    }
public:
    static void once() {
        static RegisterPlutoResources instance;
    }
};
}  // namespace

void init() {
    [[maybe_unused]] static bool initialized = []() {
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


std::string& thread_local_label() {
    thread_local static std::string label_;
    return label_;
}

std::string_view get_label() {
    return thread_local_label();
}

void set_label(std::string_view s) {
    thread_local_label().assign(s.data(),s.size());
}

void unset_label() {
    thread_local_label().clear();
}


}  // namespace pluto

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
    RegisterPlutoResources::once();
    return &MemoryResourceRegistry::instance().get(name);
}

std::string_view get_registered_name(void* mr) {
    RegisterPlutoResources::once();
    return MemoryResourceRegistry::instance().name(mr);
}

bool has_registered_resource(std::string_view name) {
    RegisterPlutoResources::once();
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

}  // namespace pluto
