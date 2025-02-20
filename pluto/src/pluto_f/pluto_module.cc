#include <string_view>
#include <cstddef>

#include "pluto/pluto.h"

extern "C" {
void c_pluto_host_set_default_resource(const char* name, int name_size) {
    pluto::host::set_default_resource( std::string_view{name,static_cast<std::size_t>(name_size)} );
}
void c_pluto_device_set_default_resource(const char* name, int name_size) {
    pluto::device::set_default_resource( std::string_view{name,static_cast<std::size_t>(name_size)} );
}
void c_pluto_scope_push() {
    pluto::scope::push();
}
void c_pluto_scope_pop() {
    pluto::scope::pop();
}
pluto::memory_resource* c_pluto_host_get_default_resource() {
    return pluto::host::get_default_resource();
}
pluto::memory_resource* c_pluto_device_get_default_resource() {
    return pluto::device::get_default_resource();
}
void* c_pluto_memory_resource_allocate(pluto::memory_resource* memory_resource, std::size_t bytes, std::size_t alignment) {
    if (alignment) {
        return memory_resource->allocate(bytes, alignment);
    }
    else {
        return memory_resource->allocate(bytes);
    }
}
void c_pluto_memory_resource_deallocate(pluto::memory_resource* memory_resource, void* memory, std::size_t bytes, std::size_t alignment) {
    if (alignment) {
        memory_resource->deallocate(memory, bytes, alignment);
    }
    else {
        memory_resource->deallocate(memory, bytes);
    }
}

}
