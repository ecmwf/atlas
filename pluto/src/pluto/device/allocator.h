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


#include "hic/hic.h"

#include "pluto/memory_resource.h"
#include "pluto/wait.h"
#include "pluto/stream.h"

namespace pluto::device {

// --------------------------------------------------------------------------------------------------------


    template <class U, class... Args>
    HIC_GLOBAL void new_on_device(U* p, Args... args) {
        // printf("new_on_device %f\n",args...);
        new (p) U(args...);
    }


    template <class U>
    HIC_GLOBAL void delete_on_device(U* p) {
        p->~U();
    }

template<typename T>
class allocator : public pluto::allocator<T> {
public:
    using value_type = T;

    allocator(memory_resource* mr, const stream& s) :
        pluto::allocator<T>::allocator(mr),
        stream_(s) {}

    allocator() :
        allocator(get_default_resource(), get_current_stream()) {}
    
    allocator(const allocator& other) :
        allocator(other.resource(), other.stream_) {}

    allocator(memory_resource* mr) :
        allocator(mr, get_current_stream()) {}

    allocator(const stream& s) :
        allocator(get_default_resource(), s) {}

    value_type* allocate(std::size_t size) {
        return pluto::allocator<T>::allocate_async(size, stream_);
    }

    void deallocate(value_type* ptr, std::size_t size) {
        pluto::allocator<T>::deallocate_async(ptr, size, stream_);
    }

    template <class U, class... Args>
    void construct(U* p, Args&&... args) {
#if HIC_COMPILER
        new_on_device<<<1, 1, 0, stream_.value<hicStream_t>()>>>(p, std::forward<Args>(args)...);
        pluto::wait(stream_);
#else
        new_on_device(p, args...);
#endif
    }

    template <class U>
    void destroy(U* p) {
#if HIC_COMPILER
        delete_on_device<<<1, 1, 0, stream_.value<hicStream_t>()>>>(p);
        pluto::wait(stream_);
#else
        delete_on_device(p);
#endif
    }
private:
    const stream& stream_;
};

// --------------------------------------------------------------------------------------------------------

}

