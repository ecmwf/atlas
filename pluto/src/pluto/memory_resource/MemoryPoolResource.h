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
#include <mutex>
#include <vector>

#include "pluto/memory_resource.h"

#include "pluto/trace.h"

namespace pluto {

pool_options get_default_pool_options();
void set_default_pool_options(pool_options);

class MemoryPoolResource : public memory_pool_resource {
public:
    MemoryPoolResource(const pool_options& options, memory_resource* upstream):
        options_(options), upstream_(upstream) {}
    MemoryPoolResource(const pool_options& options, std::unique_ptr<memory_resource>&& upstream):
        options_(options), owned_upstream_(std::move(upstream)), upstream_(owned_upstream_.get()) {}
    MemoryPoolResource(memory_resource* upstream): MemoryPoolResource(get_default_pool_options(), upstream) {}
    MemoryPoolResource(std::unique_ptr<memory_resource>&& upstream):
        MemoryPoolResource(get_default_pool_options(), std::move(upstream)) {}
    MemoryPoolResource(const pool_options& options): MemoryPoolResource(options, get_default_resource()) {}
    MemoryPoolResource(): MemoryPoolResource(get_default_pool_options(), get_default_resource()) {}

    MemoryPoolResource(memory_resource* upstream, const std::string& name, memory_tracker* memory_tracker): MemoryPoolResource(get_default_pool_options(), upstream) {
        name_ = name;
        memory_tracker_ = memory_tracker;
    }

    virtual ~MemoryPoolResource() { do_release(); }

    void release() override;

    std::size_t size() const override;

    std::size_t capacity() const override;

    void reserve(std::size_t bytes) override;

    memory_resource* upstream_resource() const override { return upstream_; }

    pool_options options() const override { return options_; }

    void do_deallocate_(void* ptr, std::size_t bytes, std::size_t alignment, bool in_callback = false);

protected:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override;
    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override;
    void* do_allocate_async(std::size_t bytes, std::size_t alignment, stream_view) override;
    void do_deallocate_async(void* ptr, std::size_t bytes, std::size_t alignment, stream_view) override;
    bool do_is_equal(const memory_resource& other) const noexcept override;
    void do_release();
    friend void callback_deallocate_async(void* stream);

    // A suitable pool or upstream resource to allocate given bytes
    memory_resource* resource(std::size_t bytes);

private:
    pool_options options_;
    std::unique_ptr<memory_resource> owned_upstream_;
    memory_resource* upstream_;
    std::string name_;
    memory_tracker* memory_tracker_{nullptr};
    std::vector<std::unique_ptr<memory_resource>> pools_;
    std::vector<std::size_t> pool_block_sizes_;
    memory_resource* pool_;
    std::size_t pool_block_size_;
    mutable std::mutex mtx_;
};

}  // namespace pluto
