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

// ------------------------------------------------------------------------------------

// These are missing includes and dummy types needed for Gator
#include <cstddef>
#include <functional>
#include <string>
#include <vector>
#include <iostream>
#include <mutex>
#include <list>

namespace pluto::yakl {
    template <typename T1, typename T2>
    inline constexpr void verbose_inform(T1, T2) {}
    struct Event{
        bool completed() const { return true; }
        bool operator ==(const Event& ) const { return true; }
    };
    inline constexpr bool yakl_mainproc() { return false; }
}

#define __YAKL_NAMESPACE_WRAPPER_BEGIN__ namespace pluto {
#define __YAKL_NAMESPACE_WRAPPER_END__ }
#include "pluto/memory_resource/detail/yakl/YAKL_Gator.h"

// ------------------------------------------------------------------------------------

#include <functional>
#include <cstddef>
#include <memory>
#include <cmath>
#include <limits>

#include "pluto/memory_resource.h"
#include "pluto/alignment.h"
#include "pluto/stream.h"

namespace pluto {

struct GatorOptions {
    std::size_t initial_size = 0;
    std::size_t grow_size    = 0;
    std::size_t block_size   = 0;
};

class GatorMemoryResource : public memory_resource {
public:

	GatorMemoryResource() : 
		GatorMemoryResource(get_default_resource()) {}

	GatorMemoryResource(memory_resource* upstream) :
		upstream_(upstream) {
		init(GatorOptions());
	}

	GatorMemoryResource(std::unique_ptr<memory_resource>&& upstream) :
		owned_upstream_(std::move(upstream)),
		upstream_(owned_upstream_.get()) {
		init(GatorOptions());
	}

	GatorMemoryResource(const GatorOptions& options) : 
		GatorMemoryResource(options, get_default_resource()) {}

	GatorMemoryResource(const GatorOptions& options, memory_resource* upstream) :
		upstream_(upstream) {
		init(options);
	}

	GatorMemoryResource(const GatorOptions& options, std::unique_ptr<memory_resource>&& upstream) :
		owned_upstream_(std::move(upstream)),
		upstream_(owned_upstream_.get()) {
		init(options);
	}

	~GatorMemoryResource() {
		if (gator_) {
			bytes_ = std::numeric_limits<std::size_t>::max();
			gator_->finalize();
		}
	}

	yakl::Gator& gator() {
		if (not gator_) {
			gator_ = std::make_unique<yakl::Gator>();
			gator_->init(
				gator_allocate_,
				gator_deallocate_,
				gator_zero_,
				initial_size_,
				grow_size_,
				block_size_
			);
		}
		return *gator_;
	}

protected:

	void init(const GatorOptions& options) {
        gator_allocate_   = [this](std::size_t bytes) { return upstream_->allocate(bytes, alignment_); };
		gator_deallocate_ = [this](void* ptr) { return upstream_->deallocate(ptr, bytes_, alignment_); };
		gator_zero_       = [](void *ptr, std::size_t bytes){};

		// Allocation alignment is guaranteed when alignment_ is Gator's blockSize
		constexpr std::size_t Gb = 1024*1024*1024;
		block_size_   = (options.block_size   ? options.block_size   : default_alignment()); //Gator default: 16*sizeof(std::size_t) = 64;
		grow_size_    = (options.grow_size    ? options.grow_size    : 1*Gb);
		initial_size_ = (options.initial_size ? options.initial_size : 1*Gb);
		alignment_    = block_size_;
	}

    void* do_allocate(std::size_t bytes, std::size_t alignment) override {
		if (std::remainder(alignment_,alignment)!=0) {
			throw std::runtime_error("alignment " + std::to_string(alignment) + " not supported. Gator pool block_size ("
				+ std::to_string( block_size_ ) + ") must multiple of requested alignment");
		}
		return gator().allocate(bytes);
    }

    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override {
		bytes_ = bytes;
		gator().free(ptr);
    }
 
    bool do_is_equal(const memory_resource& other) const noexcept override {
        if (this == &other) { return true; }
        return false;
    }

private:
	std::unique_ptr<memory_resource> owned_upstream_;
	memory_resource* upstream_;
	std::size_t bytes_; // needs to be bound to allocate/deallocate
	std::size_t alignment_; // needs to be bound to allocate/deallocate
	std::size_t initial_size_;
	std::size_t grow_size_;
	std::size_t block_size_;
	std::unique_ptr<yakl::Gator> gator_;
	std::function<void* (std::size_t)> gator_allocate_;
    std::function<void (void*)> gator_deallocate_;
	std::function<void (void*, std::size_t)> gator_zero_;
};

}
