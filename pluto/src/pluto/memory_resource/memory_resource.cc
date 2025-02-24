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
void init() {
    static std::once_flag flag;
    std::call_once(flag, [](){
#if PLUTO_DEBUGGING
        std::cout << "pluto::init()" << std::endl;
#endif
        register_resource("pluto::null_memory_resource", null_memory_resource());
        register_resource("pluto::new_delete_resource", new_delete_resource());
        register_resource("pluto::device_resource", device_resource());
        register_resource("pluto::managed_resource", managed_resource());
        register_resource("pluto::pinned_resource", pinned_resource());
        register_resource("pluto::pool_resource", pool_resource());
        register_resource("pluto::pinned_pool_resource", pinned_pool_resource());
        register_resource("pluto::device_pool_resource", device_pool_resource());
        register_resource("pluto::managed_pool_resource", managed_pool_resource());

        if (char* env = std::getenv("PLUTO_HOST_MEMORY_RESOURCE")) {
            pluto::host::set_default_resource(env);
        }
        if (char* env = std::getenv("PLUTO_DEVICE_MEMORY_RESOURCE")) {
            pluto::device::set_default_resource(env);
        }
    });
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

// class memory_resource_adaptor : public memory_resource {
// public:

//     memory_resource_adaptor(STD_PMR::memory_resource* mr, properties_t&& p): 
//         mr_(mr){
//             //properties = std::move(p);
//         }

//     STD_PMR::memory_resource* upstream_resource() { return mr_; }

// protected:
//     void* do_allocate(std::size_t bytes, std::size_t alignment) override {
//         return mr_->allocate(bytes, alignment);
//     }

//     void do_deallocate(void* p, std::size_t bytes, std::size_t alignment) override {
//         mr_->deallocate(p, bytes, alignment);
//     }

//     bool do_is_equal(const memory_resource_base& other) const noexcept override {
//         return mr_->is_equal(other);
//     }

// private:
//     friend memory_resource* get_default_resource();
//     STD_PMR::memory_resource* mr_{nullptr};
// };

memory_resource* null_memory_resource() {
    return STD_PMR::null_memory_resource();
    // static memory_resource_adaptor mr(STD_PMR::null_memory_resource(), memory_resource::properties_t{
    //     .name="pluto::null_memory_resource",
    //     .is_host_accessible=false,
    //     .is_device_accessible=false});
    // return &mr;
}
memory_resource* new_delete_resource() {
    return STD_PMR::new_delete_resource();
    // static memory_resource_adaptor mr(STD_PMR::new_delete_resource(), memory_resource::properties_t{
    //     .name="pluto::new_delete_resource",
    //     .is_host_accessible=true,
    //     .is_device_accessible=false});
    // return &mr;
}
memory_resource* get_default_resource() {
    init();
    return STD_PMR::get_default_resource();

    // auto* mr_default = dynamic_cast<memory_resource*>(std_pmr_default);
    // if (mr_default) {
    //     return mr_default;
    // }
    // if( dynamic_cast<STD_PMR::new_delete_resource_t*>(std_pmr_default) ) {
    //     return new_delete_resource();
    // }
    // if( dynamic_cast<STD_PMR::null_memory_resource_t*>(std_pmr_default) ) {
    //     return null_memory_resource();
    // }
    // static memory_resource_adaptor mr(std_pmr_default, memory_resource::properties_t{
    //     .name="std::pmr::memory_resource",
    //     .is_host_accessible=true,
    //     .is_device_accessible=false});
    // mr.mr_ = std_pmr_default;
    // return &mr;
}
void set_default_resource(STD_PMR::memory_resource* mr) {
    init();
    STD_PMR::set_default_resource(mr);
}

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
    init();
    return &MemoryResourceRegistry::instance().get(name);
}

std::string_view get_registered_name(void* mr) {
    init();
    return MemoryResourceRegistry::instance().name(mr);
}

bool has_registered_resource(std::string_view name) {
    init();
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
