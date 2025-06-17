#include <cstddef>
#include <string_view>

#include "pluto/pluto.h"

extern "C" {
void c_pluto_host_set_default_resource_name(const char* name, int name_size) {
    pluto::host::set_default_resource(std::string_view{name, static_cast<std::size_t>(name_size)});
}
void c_pluto_device_set_default_resource_name(const char* name, int name_size) {
    pluto::device::set_default_resource(std::string_view{name, static_cast<std::size_t>(name_size)});
}
void c_pluto_host_set_default_resource_ptr(pluto::memory_resource* memory_resource) {
    pluto::host::set_default_resource(memory_resource);
}
void c_pluto_device_set_default_resource_ptr(pluto::memory_resource* memory_resource) {
    pluto::device::set_default_resource(memory_resource);
}
void c_pluto_scope_push() {
    pluto::scope::push();
}
void c_pluto_scope_pop() {
    pluto::scope::pop();
}
void c_pluto_trace_enable(int enable) {
    pluto::trace::enable(enable);
}
int c_pluto_trace_enabled() {
    return pluto::trace::enabled();
}

pluto::memory_resource* c_pluto_host_get_default_resource() {
    return pluto::host::get_default_resource();
}
pluto::memory_resource* c_pluto_device_get_default_resource() {
    return pluto::device::get_default_resource();
}
void* c_pluto_memory_resource_allocate(pluto::memory_resource* memory_resource, std::size_t bytes,
                                       std::size_t alignment) {
    if (alignment) {
        return memory_resource->allocate(bytes, alignment);
    }
    else {
        return memory_resource->allocate(bytes);
    }
}
void c_pluto_memory_resource_deallocate(pluto::memory_resource* memory_resource, void* memory, std::size_t bytes,
                                        std::size_t alignment) {
    if (alignment) {
        memory_resource->deallocate(memory, bytes, alignment);
    }
    else {
        memory_resource->deallocate(memory, bytes);
    }
}
std::size_t c_pluto_memory_pool_resource_size(const pluto::memory_resource* memory_resource) {
    if (auto* pool = dynamic_cast<const pluto::memory_pool_resource*>(memory_resource)) {
        return pool->size();
    }
    return 0;
}
std::size_t c_pluto_memory_pool_resource_capacity(const pluto::memory_resource* memory_resource) {
    if (auto* pool = dynamic_cast<const pluto::memory_pool_resource*>(memory_resource)) {
        return pool->capacity();
    }
    return 0;
}
void c_pluto_memory_pool_resource_release(pluto::memory_resource* memory_resource) {
    if (auto* pool = dynamic_cast<pluto::memory_pool_resource*>(memory_resource)) {
        return pool->release();
    }
}
void c_pluto_memory_pool_resource_reserve(pluto::memory_resource* memory_resource, std::size_t bytes) {
    if (auto* pool = dynamic_cast<pluto::memory_pool_resource*>(memory_resource)) {
        return pool->reserve(bytes);
    }
}

pluto::memory_resource* c_pluto_get_registered_resource(const char* name, int name_size) {
    return pluto::get_registered_resource(std::string_view{name, static_cast<std::size_t>(name_size)});
}
pluto::memory_resource* c_pluto_new_delete_resource() {
    return pluto::new_delete_resource();
}
pluto::memory_resource* c_pluto_null_memory_resource() {
    return pluto::null_memory_resource();
}
pluto::memory_resource* c_pluto_host_resource() {
    return pluto::host_resource();
}
pluto::memory_resource* c_pluto_pinned_resource() {
    return pluto::pinned_resource();
}
pluto::memory_resource* c_pluto_device_resource() {
    return pluto::device_resource();
}
pluto::memory_resource* c_pluto_managed_resource() {
    return pluto::managed_resource();
}
pluto::memory_pool_resource* c_pluto_host_pool_resource() {
    return pluto::host_pool_resource();
}
pluto::memory_pool_resource* c_pluto_pinned_pool_resource() {
    return pluto::pinned_pool_resource();
}
pluto::memory_pool_resource* c_pluto_device_pool_resource() {
    return pluto::device_pool_resource();
}
pluto::memory_pool_resource* c_pluto_managed_pool_resource() {
    return pluto::managed_pool_resource();
}
int c_pluto_devices() {
    return pluto::devices();
}

}
