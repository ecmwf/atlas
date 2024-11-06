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

#include "pluto/memory_resource/memory_resource.h"
#include "pluto/offload/wait.h"
#include "pluto/offload/Stream.h"

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
class allocator {
public:
    using value_type = T;

    allocator(memory_resource* mr, const Stream& stream) :
        memory_resource_(mr),
        stream_(stream) {}

    allocator() :
        allocator(get_default_resource(), get_default_stream()) {}
    
    allocator(const allocator& other) :
        allocator(other.memory_resource_, other.stream_) {}

    allocator(memory_resource* mr) :
        allocator(mr, get_default_stream()) {}

    allocator(const Stream& stream) :
        allocator(get_default_resource(), stream) {}

    value_type* allocate(std::size_t size) {
        DefaultStream scope{stream_};
        return static_cast<value_type*>(memory_resource_->allocate(size * sizeof(value_type), 256));
    }

    void deallocate(value_type* ptr, std::size_t size) {
        DefaultStream scope{stream_};
        memory_resource_->deallocate(ptr, size * sizeof(value_type), 256);
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
    memory_resource* memory_resource_{nullptr};
    const Stream& stream_;
};

// --------------------------------------------------------------------------------------------------------

}

