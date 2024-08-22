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
#include <memory_resource>
#include <string>
#include <string_view>

namespace pluto {

using memory_resource = std::pmr::memory_resource;
class memory_pool_resource : public memory_resource {
public:
	virtual std::size_t size() const = 0;
	virtual std::size_t capacity() const = 0;
    virtual void release() = 0;
	virtual std::pmr::memory_resource* upstream_resource() const = 0;
    virtual std::pmr::pool_options options() const = 0;
};

inline memory_resource* null_memory_resource() { return std::pmr::null_memory_resource(); }
inline memory_resource* new_delete_resource()  { return std::pmr::new_delete_resource();  }
memory_pool_resource* pool_resource();

template<typename T>
using allocator = std::pmr::polymorphic_allocator<T>;

void init();

// --------------------------------------------------------------------------------------------------------

std::pmr::memory_resource* register_resource( std::string_view name, std::pmr::memory_resource* mr );

std::pmr::memory_resource* register_resource( std::string_view name, std::unique_ptr<std::pmr::memory_resource>&& mr );

void unregister_resource( std::string_view name );

void unregister_resources();

std::pmr::memory_resource* get_registered_resource( std::string_view name );

bool has_registered_resource(std::string_view name);

// --------------------------------------------------------------------------------------------------------

class [[nodiscard]] Register {
public:
    Register(std::string_view name, std::pmr::memory_resource* memory_resource) :
        name_(name),
        memory_resource_(register_resource(name_, memory_resource)) {}

    Register(std::string_view name, std::unique_ptr<std::pmr::memory_resource>&& memory_resource) :
        name_(name),
        memory_resource_(register_resource(name_, std::move(memory_resource))) {}

    ~Register() {
        unregister_resource(name_);
    }
    std::pmr::memory_resource* memory_resource() { return memory_resource_; }
    operator std::pmr::memory_resource*() { return memory_resource_; }

private:
    std::string name_;
    std::pmr::memory_resource* memory_resource_;
};

// --------------------------------------------------------------------------------------------------------

}
