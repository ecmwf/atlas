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
#include <cstdlib>

namespace pluto::compat {

class memory_resource {
public:
    virtual ~memory_resource() = default;
    void* allocate(std::size_t bytes, std::size_t alignment = alignof(std::max_align_t)) {
        return do_allocate(bytes, alignment);
    }
    void deallocate(void* ptr, std::size_t bytes, std::size_t alignment = alignof(std::max_align_t)) {
        do_deallocate(ptr, bytes, alignment);
    }
    bool is_equal(const memory_resource& other) const{
        return do_is_equal(other);
    }
private:
    virtual void* do_allocate(std::size_t bytes, std::size_t alignment) = 0;
    virtual void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) = 0;
    virtual bool do_is_equal(const memory_resource& other) const = 0;
};
inline bool operator==(const memory_resource& a,
                       const memory_resource& b ) noexcept {
    return &a == &b || a.is_equal(b);
}
inline bool operator!=(const memory_resource& a,
                       const memory_resource& b ) noexcept {
    return !(a == b);
}


memory_resource* get_default_resource();
memory_resource* set_default_resource(memory_resource* mr);
memory_resource* null_memory_resource();
memory_resource* new_delete_resource();

struct pool_options {
    std::size_t max_blocks_per_chunk = 0;
    std::size_t largest_required_pool_block = 0;
};

template <typename T>
class polymorphic_allocator {
public:
    using value_type = T;
    polymorphic_allocator() : mr_(get_default_resource()) {}
    polymorphic_allocator(memory_resource* mr ) : mr_(mr) {}
    polymorphic_allocator(const polymorphic_allocator& other) = default;
    template <class U>
    polymorphic_allocator(const polymorphic_allocator<U>& other) noexcept :
        mr_(other.mr_) {
    }

    value_type* allocate(std::size_t size) {
        return (value_type*)mr_->allocate(size * sizeof(value_type), alignof(value_type));
    }
    void deallocate(value_type* p, std::size_t size) {
        mr_->deallocate(p, size * sizeof(value_type), alignof(value_type));
    }
    template< class U, class... Args >
    void construct( U* p, Args&&... args ) {
        ::new(p) U(args...);
    }
    template<class U>
    void destroy(U* p) {
        p->~U();
    }
    memory_resource* resource() const {
        return mr_;
    }
private:
    memory_resource* mr_;
};


}
