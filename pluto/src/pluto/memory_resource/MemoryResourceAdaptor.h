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

#include <functional>

namespace pluto {

// --------------------------------------------------------------------------------------------------------

class MemoryResourceAdaptor : public memory_resource {
public:
    template <typename Allocate, typename Deallocate>
    MemoryResourceAdaptor(Allocate allocate, Deallocate deallocate):
        callback_allocate_(allocate), callback_deallocate_(deallocate) {}

protected:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override {
        return callback_allocate_(bytes, alignment);
    }

    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override {
        return callback_deallocate_(ptr, bytes, alignment);
    }

    bool do_is_equal(const memory_resource& other) const noexcept override {
        if (this == &other) {
            return true;
        }
        return false;
    }

private:
    std::function<void*(std::size_t, std::size_t)> callback_allocate_;
    std::function<void(void*, std::size_t, std::size_t)> callback_deallocate_;
};

// --------------------------------------------------------------------------------------------------------

}  // namespace pluto
