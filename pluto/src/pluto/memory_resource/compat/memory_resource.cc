/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "pluto/memory_resource/compat/memory_resource"

#include <cstdlib> // posix_memalign
#include <exception>

#include <iostream>

namespace pluto::compat {

class new_delete_resource_t : public memory_resource {
public:
    new_delete_resource_t() {}
protected:
    virtual void* do_allocate(std::size_t bytes, std::size_t alignment) override {

        // This implementation is backward compatible, not relying on C++17 headers,
        // which is the point of this compat in the first place
        void* ptr;
        alignment = std::max(alignment, alignof(std::max_align_t));
        int err = ::posix_memalign((void**)&ptr, alignment, bytes);
        if (err) {
            throw std::bad_alloc();
        }
        return ptr;

        // This implementation requires C++17 header
        // void* ptr = ::operator new(bytes, std::align_val_t(alignment));
        // return ptr;
    }
    virtual void do_deallocate(void* ptr, std::size_t /*bytes*/, std::size_t /*alignment*/) override {

        // This implementation is backward compatible, not relying on C++17 headers,
        // which is the point of this compat in the first place
        ::free(ptr);

        // This implementation requires C++17 header
        // ::operator delete(ptr, bytes, std::align_val_t(alignment));
    }
    virtual bool do_is_equal(const memory_resource& other) const noexcept override {
        return &other == this;
    }
};

class null_memory_resource_t : public memory_resource {
public:
    null_memory_resource_t() {}
protected:
    virtual void* do_allocate(std::size_t /*bytes*/, std::size_t /*alignment*/) override {
        throw std::bad_alloc();
    }
    virtual void do_deallocate(void* /*ptr*/, std::size_t /*bytes*/, std::size_t /*alignment*/) override {}
    virtual bool do_is_equal(const memory_resource& other) const noexcept override {
        return &other == this;
    }
};


template<typename T>
struct constant_init {
    union {
        T obj;
    };
    constexpr constant_init() : obj() { }

    template<typename U>
    explicit constexpr constant_init(U arg) : obj(arg) { }

    ~constant_init() { /* do nothing, union member is not destroyed */ }
};

constant_init<new_delete_resource_t>  newdel_res{};
constant_init<null_memory_resource_t> null_res{};

class default_resource {
public:
    static memory_resource* get() {
        return instance().mr_;
    }
    static memory_resource* set(memory_resource* mr) {
        auto previous = instance().mr_;
        instance().mr_ = mr;
        return previous;
    }
private:
    static default_resource& instance() {
        static default_resource instance;
        return instance;
    }
    default_resource () {
        mr_ = new_delete_resource();
    }
    memory_resource* mr_;
};

memory_resource* get_default_resource() {
    return default_resource::get();
}
memory_resource* set_default_resource(memory_resource* mr) {
    return default_resource::set(mr);
}
memory_resource* null_memory_resource() {
    return &null_res.obj;
}
memory_resource* new_delete_resource() {
    return &newdel_res.obj;
}


}
