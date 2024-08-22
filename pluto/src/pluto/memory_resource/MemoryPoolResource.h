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
#include <memory_resource>
#include <vector>

#include "pluto/memory_resource/memory_resource.h"

namespace pluto {

std::pmr::pool_options get_default_pool_options();
void set_default_pool_options(std::pmr::pool_options);

class MemoryPoolResource : public memory_pool_resource {
public:
	MemoryPoolResource(const std::pmr::pool_options& options, std::pmr::memory_resource* upstream) :
		options_(options),
		upstream_(upstream) {
	}
	MemoryPoolResource(const std::pmr::pool_options& options, std::unique_ptr<std::pmr::memory_resource>&& upstream) :
		options_(options),
		owned_upstream_(std::move(upstream)),
		upstream_(owned_upstream_.get()) {
	}
	MemoryPoolResource(std::pmr::memory_resource* upstream) :
		MemoryPoolResource(get_default_pool_options(),upstream) {
	}
	MemoryPoolResource(std::unique_ptr<std::pmr::memory_resource>&& upstream) :
		MemoryPoolResource(get_default_pool_options(), std::move(upstream)) {
	}
	MemoryPoolResource(const std::pmr::pool_options& options) :
		MemoryPoolResource(options, std::pmr::get_default_resource()) {
	}
	MemoryPoolResource() :
		MemoryPoolResource(get_default_pool_options(), std::pmr::get_default_resource()) {
	}

	void release() override {
		pools_.clear();
	}

	std::size_t size() const override;

	std::size_t capacity() const override;

	std::pmr::memory_resource* upstream_resource() const override {
		return upstream_;
	}

	std::pmr::pool_options options() const override {
		return options_;
	}

protected:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override {
		return resource(bytes)->allocate(bytes, alignment);
    }

    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override {
		resource(bytes)->deallocate(ptr, bytes, alignment);
		// possible to cleanup no longer used gators
		cleanup_unused_gators();
    }
 
    bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override {
        if (this == &other) { return true; }
        return false;
    }

	// A suitable pool or upstream resource to allocate and deallocate given bytes
	std::pmr::memory_resource* resource(std::size_t bytes);

	void cleanup_unused_gators();

private:
	std::pmr::pool_options options_;
	std::unique_ptr<std::pmr::memory_resource> owned_upstream_;
	std::pmr::memory_resource* upstream_;
	std::vector<std::unique_ptr<std::pmr::memory_resource>> pools_;
	std::vector<std::size_t> pool_block_sizes_;
	std::pmr::memory_resource* pool_;
	std::size_t pool_block_size_;
};

}
