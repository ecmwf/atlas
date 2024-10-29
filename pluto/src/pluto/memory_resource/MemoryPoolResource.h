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
#include <vector>
#include <mutex>

#include "pluto/memory_resource/memory_resource.h"

namespace pluto {

pool_options get_default_pool_options();
void set_default_pool_options(pool_options);

class MemoryPoolResource : public memory_pool_resource {
public:
	MemoryPoolResource(const pool_options& options, memory_resource* upstream) :
		options_(options),
		upstream_(upstream) {
	}
	MemoryPoolResource(const pool_options& options, std::unique_ptr<memory_resource>&& upstream) :
		options_(options),
		owned_upstream_(std::move(upstream)),
		upstream_(owned_upstream_.get()) {
	}
	MemoryPoolResource(memory_resource* upstream) :
		MemoryPoolResource(get_default_pool_options(),upstream) {
	}
	MemoryPoolResource(std::unique_ptr<memory_resource>&& upstream) :
		MemoryPoolResource(get_default_pool_options(), std::move(upstream)) {
	}
	MemoryPoolResource(const pool_options& options) :
		MemoryPoolResource(options, get_default_resource()) {
	}
	MemoryPoolResource() :
		MemoryPoolResource(get_default_pool_options(), get_default_resource()) {
	}

	void release() override {
		std::lock_guard lock(mtx_);
		pools_.clear();
	}

	std::size_t size() const override;

	std::size_t capacity() const override;

	memory_resource* upstream_resource() const override {
		return upstream_;
	}

	pool_options options() const override {
		return options_;
	}

protected:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override;
    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override;
    bool do_is_equal(const memory_resource& other) const noexcept override;

	// A suitable pool or upstream resource to allocate given bytes
	memory_resource* resource(std::size_t bytes);

private:
	pool_options options_;
	std::unique_ptr<memory_resource> owned_upstream_;
	memory_resource* upstream_;
	std::vector<std::unique_ptr<memory_resource>> pools_;
	std::vector<std::size_t> pool_block_sizes_;
	memory_resource* pool_;
	std::size_t pool_block_size_;
	mutable std::mutex mtx_;
};

}
