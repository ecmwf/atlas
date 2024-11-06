/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "MemoryPoolResource.h"

#include "GatorMemoryResource.h"
#include "TraceMemoryResource.h"
#include "pluto/util/Trace.h"

#include <cstdlib>
#include <cctype>
#include <string>
#include <map>

namespace pluto {

static std::size_t to_bytes(const std::string& str) {
    auto unit_to_bytes = [](char unit) {
       static const std::map<char, std::size_t> map{
            {'G',1024*1024*1024},
            {'M',1024*1024},
            {'K',1024},
            {'B',1}
        };
        return map.at(static_cast<char>(std::toupper(unit)));
    };
    for (char unit: {'G','g','M','m','K','k','B','b'}) {
        if (auto pos = str.find(unit); pos != std::string::npos) {
            return unit_to_bytes(unit) * std::stoull(str.substr(0,pos));
        }
    }
    return std::stoull(str);
}

static std::string bytes_to_string(std::size_t bytes) {
    std::string str;
    constexpr std::size_t MB = 1024*1024;
    constexpr std::size_t GB = 1024*MB;
    if( bytes >= GB ) {
        str = std::to_string(bytes/GB) + "GB";
    }
    else {
        str = std::to_string(bytes/MB) + "MB";
    }
    return str;
}


static pool_options default_pool_options_;
static bool default_pool_options_setup_ = false;

pool_options get_default_pool_options() {
    if (not default_pool_options_setup_) {
        if (const char* env = std::getenv("PLUTO_LARGEST_REQUIRED_POOL_BLOCK"); env != nullptr) {
            default_pool_options_.largest_required_pool_block = to_bytes(env);
        }
        if (const char* env = std::getenv("PLUTO_MAX_BLOCKS_PER_CHUNK"); env != nullptr) {
            default_pool_options_.max_blocks_per_chunk = std::atoi(env);
        }
        default_pool_options_setup_ = true;
    }
    return default_pool_options_;
}

void set_default_pool_options(pool_options options) {
    default_pool_options_ = options;
    default_pool_options_setup_ = true;
}

static GatorMemoryResource* to_gator_resource(memory_resource* pool) {
    GatorMemoryResource* gator;
    if (TraceMemoryResource* traced = dynamic_cast<TraceMemoryResource*>(pool)) {
        gator = dynamic_cast<GatorMemoryResource*>(traced->upstream_resource());
    }
    else {
        gator = dynamic_cast<GatorMemoryResource*>(pool);
    }
    return gator;
}

void* MemoryPoolResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    std::lock_guard lock(mtx_);
    return resource(bytes)->allocate(bytes, alignment);
}

void MemoryPoolResource::do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
    std::lock_guard lock(mtx_);
    if (pools_.size() == 1) {
        pools_[0]->deallocate(ptr, bytes, alignment);
    }
    else {
        for (int i=0; i<pool_block_sizes_.size(); ++i) {
            if (pools_[i]) {
                auto& gator = to_gator_resource(pools_[i].get())->gator();
                if (gator.thisIsMyPointer(ptr)) {
                    pools_[i]->deallocate(ptr, bytes, alignment);
                }

                // Cleanup empty gator when a larger gator exists
                if (gator.get_bytes_currently_allocated() == 0 && pool_block_size_ > pool_block_sizes_[i]) {
#if PLUTO_DEBUGGING
                std::cout << " - Releasing memory_pool["<<i<<"] with block_size " << bytes_to_string(pool_block_sizes_[i]) << " and capacity " << bytes_to_string(gator.get_pool_capacity()) << std::endl;
#endif
                    pools_[i].reset();
                }
            }
        }
    }
}

void MemoryPoolResource::reserve(std::size_t bytes) {
    deallocate(allocate(bytes), bytes);
}

bool MemoryPoolResource::do_is_equal(const memory_resource& other) const noexcept {
    if (this == &other) { return true; }
    return false;
}

memory_resource* MemoryPoolResource::resource(std::size_t bytes) {
    constexpr std::size_t MB = 1024*1024;
    constexpr std::size_t GB = 1024*MB;
    if (options_.largest_required_pool_block > 0 && bytes > options_.largest_required_pool_block) {
        return upstream_;
    }
    if (pools_.empty()) {
        if (options_.largest_required_pool_block != 0) {
            pool_block_sizes_ = {options_.largest_required_pool_block};
        }
        else {
            pool_block_sizes_ = {256*MB, 1*GB, 4*GB, 16*GB, 64*GB, 256*GB, 1024*GB};
        }
        pools_.resize(pool_block_sizes_.size());
        pool_block_size_ = 0;
    }

    auto upper_or_equal = std::upper_bound(pool_block_sizes_.begin(), pool_block_sizes_.end(), bytes-1);
    if (upper_or_equal == pool_block_sizes_.end()) {
        return upstream_;
    }

    if (*upper_or_equal > pool_block_size_) {
        pool_block_size_ = *upper_or_equal;
        std::size_t pool_index = upper_or_equal - pool_block_sizes_.begin();
        std::size_t blocks_per_chunk = std::max(options_.max_blocks_per_chunk, static_cast<std::size_t>(1));
        GatorOptions options;
        options.initial_size = blocks_per_chunk*pool_block_size_;
        options.grow_size    = blocks_per_chunk*pool_block_size_;
        if (TraceOptions::instance().enabled) {
            pools_[pool_index] =
                std::make_unique<TraceMemoryResource>("gator["+bytes_to_string(pool_block_size_)+"]",
                std::make_unique<GatorMemoryResource>(options, upstream_) );
        }
        else {
            pools_[pool_index] = std::make_unique<GatorMemoryResource>(options, upstream_);
        }
        pool_ = pools_[pool_index].get();
    }
    return pool_;
}

std::size_t MemoryPoolResource::size() const {
    std::lock_guard lock(mtx_);
    std::size_t _size{0};
    for (const auto& pool: pools_) {
        if (pool) {
            GatorMemoryResource* gator = to_gator_resource(pool.get());
            if (gator) {
                _size += gator->gator().get_bytes_currently_allocated();
            }
        }
    }
    return _size;
}

std::size_t MemoryPoolResource::capacity() const {
    std::lock_guard lock(mtx_);
    std::size_t _capacity{0};
    for (const auto& pool: pools_) {
        if (pool) {
            GatorMemoryResource* gator = to_gator_resource(pool.get());
            if (gator) {
                _capacity += gator->gator().get_pool_capacity();
            }
        }
    }
    return _capacity;
}


}
