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

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>

#include "pluto/pluto_config.h"

#if PLUTO_HAVE_PMR
#include <memory_resource>
#define STD_PMR std::pmr
#else
#include "pluto/memory_resource/compat/memory_resource"
#define STD_PMR pluto::compat
#endif

namespace pluto {

void init();

using memory_resource = STD_PMR::memory_resource;
using pool_options = STD_PMR::pool_options;
// template<typename T>
// using allocator = STD_PMR::polymorphic_allocator<T>;

inline memory_resource* null_memory_resource() { return STD_PMR::null_memory_resource(); }
inline memory_resource* new_delete_resource()  { return STD_PMR::new_delete_resource();  }
inline memory_resource* get_default_resource() { 
    init();
    return STD_PMR::get_default_resource(); }
inline void set_default_resource(memory_resource* mr) { STD_PMR::set_default_resource(mr); }

template<typename T>
class allocator : public STD_PMR::polymorphic_allocator<T> {
public:
    using STD_PMR::polymorphic_allocator<T>::polymorphic_allocator;
    allocator() : STD_PMR::polymorphic_allocator<T>(get_default_resource()) {}
};

class memory_pool_resource : public memory_resource {
public:
	virtual std::size_t size() const = 0;
	virtual std::size_t capacity() const = 0;
    virtual void release() = 0;
	virtual memory_resource* upstream_resource() const = 0;
    virtual pool_options options() const = 0;
    virtual void reserve(std::size_t) = 0;
};

memory_pool_resource* pool_resource();


// --------------------------------------------------------------------------------------------------------

memory_resource* register_resource( std::string_view name, memory_resource* mr );

memory_resource* register_resource( std::string_view name, std::unique_ptr<memory_resource>&& mr );

void unregister_resource( std::string_view name );

void unregister_resources();

memory_resource* get_registered_resource( std::string_view name );

bool has_registered_resource(std::string_view name);

// --------------------------------------------------------------------------------------------------------

class [[nodiscard]] Register {
public:
    Register(std::string_view name, memory_resource* mr) :
        name_(name),
        memory_resource_(register_resource(name_, mr)) {}

    Register(std::string_view name, std::unique_ptr<memory_resource>&& mr) :
        name_(name),
        memory_resource_(register_resource(name_, std::move(mr))) {}

    ~Register() {
        unregister_resource(name_);
    }
    pluto::memory_resource* memory_resource() { return memory_resource_; }
    operator pluto::memory_resource*() { return memory_resource_; }

private:
    std::string name_;
    pluto::memory_resource* memory_resource_;
};

// --------------------------------------------------------------------------------------------------------

}
