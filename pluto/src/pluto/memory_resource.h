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

#include "pluto/alignment.h"
#include "pluto/pluto_config.h"

#pragma push_macro("STD_PMR")
#undef STD_PMR

#if PLUTO_HAVE_PMR
#include <memory_resource>
#define STD_PMR std::pmr
#else
#include "pluto/memory_resource/compat/memory_resource"
#define STD_PMR pluto::compat
#endif

#include "pluto/stream.h"

namespace pluto {

std::string_view get_label();
void set_label(std::string_view);
void unset_label();
class scoped_label {
public:
    scoped_label(std::string_view s) :
        previous_(get_label()) {
        set_label(s); }
    ~scoped_label() {
        set_label(previous_);
    }
private:
    std::string previous_;
};

using memory_resource = STD_PMR::memory_resource;
using pool_options    = STD_PMR::pool_options;

inline memory_resource* null_memory_resource() {
    return STD_PMR::null_memory_resource();
}

inline memory_resource* new_delete_resource() {
    return STD_PMR::new_delete_resource();
}

inline memory_resource* get_default_resource() {
    return STD_PMR::get_default_resource();
}

inline void set_default_resource(memory_resource* mr) {
    STD_PMR::set_default_resource(mr);
}

class async_memory_resource : public memory_resource {
public:
    using memory_resource::memory_resource;

    void* allocate_async(std::size_t bytes, std::size_t alignment, stream_view s) {
        return do_allocate_async(bytes, alignment, s);
    }
    void deallocate_async(void* ptr, std::size_t bytes, std::size_t alignment, stream_view s) {
        do_deallocate_async(ptr, bytes, alignment, s);
    }

private:
    virtual void* do_allocate_async(std::size_t bytes, std::size_t alignment, stream_view)             = 0;
    virtual void do_deallocate_async(void* ptr, std::size_t bytes, std::size_t alignment, stream_view) = 0;
};

template <typename T>
class allocator : public STD_PMR::polymorphic_allocator<T> {
    using base_t = STD_PMR::polymorphic_allocator<T>;

public:
    using value_type = typename base_t::value_type;
    using base_t::polymorphic_allocator;

    allocator(): allocator(get_default_resource()) {}

    allocator(memory_resource* mr): base_t(mr), async_mr_(dynamic_cast<async_memory_resource*>(base_t::resource())) {}

    allocator(const base_t& other):
        base_t(other), async_mr_(dynamic_cast<async_memory_resource*>(base_t::resource())) {}

    template <class U>
    allocator(const STD_PMR::polymorphic_allocator<U>& other) noexcept:
        base_t(other), async_mr_(dynamic_cast<async_memory_resource*>(base_t::resource())) {}

    using base_t::allocate;
    value_type* allocate(std::string_view label, std::size_t size) {
        if (label.size()) {
            scoped_label scoped(label);
            return base_t::allocate(size);
        }
        return base_t::allocate(size);
    }

    using base_t::deallocate;
    void deallocate(std::string_view label, value_type* ptr, std::size_t size) {
        if (label.size()) {
            scoped_label scoped(label);
            base_t::deallocate(ptr, size);
        }
        else {
            base_t::deallocate(ptr, size);
        }
    }

    value_type* allocate_async(std::size_t size, stream_view s) {
        if (async_mr_) {
            return (value_type*)async_mr_->allocate_async(size * sizeof(value_type), pluto::default_alignment(), s);
        }
        else {
            scoped_stream scoped(s);
            return base_t::allocate(size);
        }
    }

    void deallocate_async(value_type* ptr, std::size_t size, stream_view s) {
        if (async_mr_) {
            async_mr_->deallocate_async(ptr, size * sizeof(value_type), pluto::default_alignment(), s);
        }
        else {
            scoped_stream scoped(s);
            base_t::deallocate(ptr, size);
        }
    }

    value_type* allocate_async(std::string_view label, std::size_t size, stream_view s) {
        if (label.empty()) {
            return allocate_async(size, s);
        }
        scoped_label scoped(label);
        return allocate_async(size, s);
    }

    void deallocate_async(std::string_view label, value_type* ptr, std::size_t size, stream_view s) {
        if (label.empty()) {
            deallocate_async(ptr, size, s);
        }
        scoped_label scoped(label);
        deallocate_async(ptr, size, s);
    }


private:
    async_memory_resource* async_mr_{nullptr};
};

class memory_pool_resource : public async_memory_resource {
public:
    virtual std::size_t size() const                   = 0;
    virtual std::size_t capacity() const               = 0;
    virtual void release()                             = 0;
    virtual memory_resource* upstream_resource() const = 0;
    virtual pool_options options() const               = 0;
    virtual void reserve(std::size_t)                  = 0;
};

// --------------------------------------------------------------------------------------------------------

memory_resource* register_resource(std::string_view name, memory_resource* mr);

memory_resource* register_resource(std::string_view name, std::unique_ptr<memory_resource>&& mr);

void unregister_resource(std::string_view name);

void unregister_resources();

memory_resource* get_registered_resource(std::string_view name);

std::string_view get_registered_name(void*);

bool has_registered_resource(std::string_view name);

// --------------------------------------------------------------------------------------------------------

class [[nodiscard]] Register {
public:
    Register(std::string_view name, memory_resource* mr): name_(name), memory_resource_(register_resource(name_, mr)) {}

    Register(std::string_view name, std::unique_ptr<memory_resource>&& mr):
        name_(name), memory_resource_(register_resource(name_, std::move(mr))) {}

    ~Register() { unregister_resource(name_); }
    pluto::memory_resource* memory_resource() { return memory_resource_; }
    operator pluto::memory_resource*() { return memory_resource_; }

private:
    std::string name_;
    pluto::memory_resource* memory_resource_;
};

// --------------------------------------------------------------------------------------------------------

}  // namespace pluto

#pragma pop_macro("STD_PMR")
