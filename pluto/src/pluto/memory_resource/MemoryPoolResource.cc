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

#include "TraceMemoryResource.h"
#include "detail/GatorMemoryResource.h"
#include "pluto/memory.h"
#include "pluto/trace.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <map>
#include <queue>
#include <string>

#include "hic/hic.h"

namespace pluto {

static std::size_t to_bytes(const std::string& str) {
    auto unit_to_bytes = [](char unit) {
        static const std::map<char, std::size_t> map{
            {'G', 1024 * 1024 * 1024}, {'M', 1024 * 1024}, {'K', 1024}, {'B', 1}};
        return map.at(static_cast<char>(std::toupper(unit)));
    };
    for (char unit : {'G', 'g', 'M', 'm', 'K', 'k', 'B', 'b'}) {
        if (auto pos = str.find(unit); pos != std::string::npos) {
            return unit_to_bytes(unit) * std::stoull(str.substr(0, pos));
        }
    }
    return std::stoull(str);
}

static std::string bytes_to_string(std::size_t bytes) {
    std::string str;
    constexpr std::size_t MB = 1024 * 1024;
    constexpr std::size_t GB = 1024 * MB;
    if (bytes >= GB) {
        str = std::to_string(bytes / GB) + "GB";
    }
    else {
        str = std::to_string(bytes / MB) + "MB";
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
    default_pool_options_       = options;
    default_pool_options_setup_ = true;
}


// static GatorMemoryResource* to_gator_resource(memory_resource* pool) {
//     memory_resource* upstream = pool;
//     GatorMemoryResource* gator;
//     do {
//         gator = dynamic_cast<GatorMemoryResource*>(upstream);
//         upstream = dynamic_cast<memory_resource*>(pool->properties.upstream_resource);
//     } while(not gator);
//     return gator;
// }

static GatorMemoryResource* to_gator_resource(memory_resource* pool) {
    GatorMemoryResource* gator;
#if PLUTO_DEBUGGING
    if (TraceMemoryResource* traced = dynamic_cast<TraceMemoryResource*>(pool)) {
        gator = dynamic_cast<GatorMemoryResource*>(traced->upstream_resource());
    }
    else {
        gator = dynamic_cast<GatorMemoryResource*>(pool);
    }
#else
    gator = dynamic_cast<GatorMemoryResource*>(pool);
#endif
    return gator;
}

void* MemoryPoolResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    std::lock_guard lock(mtx_);
    auto* mr = resource(bytes);
    void* ptr = mr->allocate(bytes, alignment);

    bool used_upstream = (*mr == *upstream_);
    if (used_upstream) {
        if (trace::enabled() && name_.size()) {
            trace::out << "PLUTO_TRACE    --> used instead of " << name_ << " as bytes > largest_required_pool_block (" << trace::format_bytes(options_.largest_required_pool_block) << ")\n";
        }
    }
    else {
        if (memory_tracker_) {
            memory_tracker_->allocate(bytes);
        }
        if (name_.size() && trace::enabled()) {
            trace::log::allocate(get_label(), ptr, bytes, alignment, name_, memory_tracker_);
        }
    }
    return ptr;
}

void* MemoryPoolResource::do_allocate_async(std::size_t bytes, std::size_t alignment, stream_view s) {
    std::lock_guard lock(mtx_);
    auto* mr       = resource(bytes);
    auto* async_mr = dynamic_cast<async_memory_resource*>(mr);
    void* ptr;
    if (async_mr) {
        ptr = async_mr->allocate_async(bytes, alignment, s);
    }
    else {
        ptr = mr->allocate(bytes, alignment);
    }
    bool used_upstream = (*mr == *upstream_);
    if (not used_upstream) {
        if (memory_tracker_) {
            memory_tracker_->allocate(bytes);
        }
        if (name_.size() && trace::enabled()) {
            trace::log::allocate_async(get_label(), ptr, bytes, alignment, s.value(), name_, memory_tracker_);
        }
    }
    return ptr;
}

void MemoryPoolResource::do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
    do_deallocate_(ptr, bytes, alignment);
}

void MemoryPoolResource::do_deallocate_(void* ptr, std::size_t bytes, std::size_t alignment, bool in_callback) {
    std::lock_guard lock(mtx_);

    bool use_upstream = false;
    if (options_.largest_required_pool_block > 0 && bytes > options_.largest_required_pool_block) {
        use_upstream = true;
    }

    if (not in_callback && not use_upstream) {
        if (memory_tracker_) {
            memory_tracker_->deallocate(bytes);
        }
        if (trace::enabled() && name_.size()) {
            trace::log::deallocate(get_label(), ptr, bytes, alignment, name_, memory_tracker_);
        }
    }

    if (use_upstream) {
        upstream_->deallocate(ptr, bytes, alignment);
        if (trace::enabled() && name_.size()) {
            trace::out << "PLUTO_TRACE    --> used instead of " << name_ << " as bytes > largest_required_pool_block (" << trace::format_bytes(options_.largest_required_pool_block) << ")\n";
        }
    }
    else if (pools_.size() == 1) {
        pools_[0]->deallocate(ptr, bytes, alignment);
    }
    else {
        for (std::size_t i = 0; i < pool_block_sizes_.size(); ++i) {
            if (pools_[i]) {
                auto& gator = to_gator_resource(pools_[i].get())->gator();
                if (gator.thisIsMyPointer(ptr)) {
                    pools_[i]->deallocate(ptr, bytes, alignment);
                }

                // Cleanup empty gator when a larger gator exists
                if (gator.get_bytes_currently_allocated() == 0 && pool_block_size_ > pool_block_sizes_[i]) {
#if PLUTO_DEBUGGING
                    std::cout << " - Releasing memory_pool[" << i << "] with block_size "
                              << bytes_to_string(pool_block_sizes_[i]) << " and capacity "
                              << bytes_to_string(gator.get_pool_capacity()) << std::endl;
#endif
                    pools_[i].reset();
                }
            }
        }
    }
}


struct AsyncAllocData {
    MemoryPoolResource* resource;
    void* ptr;
    std::size_t bytes;
    std::size_t alignment;
};

std::map<void*, std::queue<AsyncAllocData>> stream_callback_queue_;

void callback_deallocate_async(void* stream) {
    auto& stream_queue        = stream_callback_queue_[stream];
    const AsyncAllocData& ctx = stream_queue.front();
    constexpr bool in_callback = true;
    ctx.resource->do_deallocate_(ctx.ptr, ctx.bytes, ctx.alignment, in_callback);
    stream_queue.pop();
}

void MemoryPoolResource::do_deallocate_async(void* ptr, std::size_t bytes, std::size_t alignment, stream_view s) {
    // stream.wait();
    // Wait for stream to finish for safety.
    // We should not deallocate data when it may still be in use in the stream!
    // TODO: implement using
    //    __host__ ​cudaError_t cudaLaunchHostFunc ( cudaStream_t stream, cudaHostFn_t fn, void* userData )
    auto& stream_queue = stream_callback_queue_[s.value()];
    stream_queue.emplace(AsyncAllocData{this, ptr, bytes, alignment});
    HIC_CALL(hicLaunchHostFunc(s.value<hicStream_t>(), callback_deallocate_async, s.value()));

    if (memory_tracker_) {
        memory_tracker_->deallocate(bytes);
    }
    if (name_.size() && trace::enabled()) {
        trace::log::deallocate_async(get_label(),ptr, bytes, alignment, s.value(), name_, memory_tracker_);
    }
}


void MemoryPoolResource::reserve(std::size_t bytes) {
    // This reserve uses the pool api to allocate a large chunk of bytes
    // This could fail when the largest_required_pool_block is set to a lower value than bytes

    // Disable memory tracking and tracing for the pool itself, it should not count as a pool allocation.
    // Note that the upstream memory resource may still be tracking!
    bool memory_tracker_previously_enabled = memory_tracker_ ? memory_tracker_->enable(false) : false ;
    std::string stored_name{name_};
    name_.clear(); // a trick to disable logging this as an allocation of a variable

    deallocate(allocate(bytes), bytes);

    // Restore memory tracking and tracing
    name_ = stored_name; // restore name so that logging of variables will be available again
    if (memory_tracker_previously_enabled) {
        memory_tracker_->enable(true);
    }
}

bool MemoryPoolResource::do_is_equal(const memory_resource& other) const noexcept {
    return (this == &other);
}

memory_resource* MemoryPoolResource::resource(std::size_t bytes) {
    constexpr std::size_t MB = 1024 * 1024;
    constexpr std::size_t GB = 1024 * MB;
    if (options_.largest_required_pool_block > 0 && bytes > options_.largest_required_pool_block) {
        return upstream_;
    }
    if (pools_.empty()) {
        if (options_.largest_required_pool_block != 0) {
            pool_block_sizes_ = {options_.largest_required_pool_block};
        }
        else {
            pool_block_sizes_ = {256 * MB, 1 * GB, 4 * GB, 16 * GB, 64 * GB, 256 * GB, 1024 * GB};
        }
        pools_.resize(pool_block_sizes_.size());
        pool_block_size_ = 0;
    }

    auto upper_or_equal = std::upper_bound(pool_block_sizes_.begin(), pool_block_sizes_.end(), bytes - 1);
    if (upper_or_equal == pool_block_sizes_.end()) {
        return upstream_;
    }

    if (*upper_or_equal > pool_block_size_) {
        pool_block_size_             = *upper_or_equal;
        std::size_t pool_index       = upper_or_equal - pool_block_sizes_.begin();
        std::size_t blocks_per_chunk = std::max(options_.max_blocks_per_chunk, static_cast<std::size_t>(1));
        GatorOptions options;
        options.initial_size = blocks_per_chunk * pool_block_size_;
        options.grow_size    = blocks_per_chunk * pool_block_size_;
        if (trace::enabled() && PLUTO_DEBUGGING) {
            pools_[pool_index] = std::make_unique<TraceMemoryResource>(
                "gator[" + bytes_to_string(pool_block_size_) + "]",
                std::make_unique<GatorMemoryResource>(options, std::make_unique<TraceMemoryResource>(upstream_)));
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
    for (const auto& pool : pools_) {
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
    for (const auto& pool : pools_) {
        if (pool) {
            GatorMemoryResource* gator = to_gator_resource(pool.get());
            if (gator) {
                _capacity += gator->gator().get_pool_capacity();
            }
        }
    }
    return _capacity;
}

void MemoryPoolResource::release() {
    std::lock_guard lock(mtx_);
#if PLUTO_DEBUGGING
    for (std::size_t i = 0; i < pool_block_sizes_.size(); ++i) {
        if (pools_[i]) {
            auto& gator = to_gator_resource(pools_[i].get())->gator();
            // Cleanup empty gator when a larger gator exists
            std::cout << " - Releasing memory_pool[" << i << "] with block_size "
                      << bytes_to_string(pool_block_sizes_[i]) << " and capacity "
                      << bytes_to_string(gator.get_pool_capacity()) << std::endl;
            pools_[i].reset();
        }
    }
#endif
    if (name_.size() && trace::enabled()) {
        trace::out << name_ << "::release()" << std::endl;
    }
    pools_.clear();
}


}  // namespace pluto
