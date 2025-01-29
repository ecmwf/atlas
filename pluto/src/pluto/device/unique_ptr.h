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

#include <memory>
#include <utility>

#include "allocator.h"

namespace pluto::device {

// --------------------------------------------------------------------------------------------------------

template <class Alloc>
class Deleter {
private:
    Alloc alloc_;

public:
    using value_type = typename Alloc::value_type;

    Deleter() = default;
    Deleter(const Alloc& alloc): alloc_(alloc) {}
    void operator()(value_type* p) {
        alloc_.destroy(p);
        alloc_.deallocate(p, 1);
    }
};

template <typename T>
using unique_ptr = std::unique_ptr<T, Deleter<allocator<T>>>;

template <class T, class... Args>
unique_ptr<T> make_unique(Args&&... args) {
    allocator<T> alloc;
    T* p = alloc.allocate(1);
    alloc.construct(p, std::forward<Args>(args)...);
    return unique_ptr<T>(p, Deleter(alloc));
}

// --------------------------------------------------------------------------------------------------------

}  // namespace pluto::device
