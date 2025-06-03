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
#include <iostream>
#include <list>
#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace pluto::yakl {
template <typename T1, typename T2>
inline constexpr void verbose_inform(T1, T2) {}
struct Event {
    bool completed() const { return true; }
    bool operator==(const Event&) const { return true; }
};
inline constexpr bool yakl_mainproc() {
    return false;
}
}  // namespace pluto::yakl

#define __YAKL_NAMESPACE_WRAPPER_BEGIN__ namespace pluto {
#define __YAKL_NAMESPACE_WRAPPER_END__ }

#if defined(__NVCOMPILER)
#define PLUTO_SUPPRESS_WARNINGS_PUSH _Pragma("diag push")
#define PLUTO_SUPPRESS_WARNINGS_POP _Pragma("diag pop")
#define PLUTO_SUPPRESS_WARNINGS_INTEGER_SIGN_CHANGE _Pragma("diag_suppress integer_sign_change")
#define PLUTO_SUPPRESS_WARNINGS_CODE_IS_UNREACHABLE _Pragma("diag_suppress code_is_unreachable")
#elif defined(__INTEL_COMPILER)
#define PLUTO_SUPPRESS_WARNINGS_PUSH _Pragma("warning push")
#define PLUTO_SUPPRESS_WARNINGS_POP _Pragma("warning pop")
#define PLUTO_SUPPRESS_WARNINGS_INTEGER_SIGN_CHANGE _Pragma("warning disable 68")
#elif defined(__GNUC__)
#define PLUTO_SUPPRESS_WARNINGS_PUSH                                               \
    _Pragma("GCC diagnostic push") _Pragma("GCC diagnostic ignored \"-Wpragmas\"") \
        _Pragma("GCC diagnostic ignored \"-Wunknown-warning-option\"")
#define PLUTO_SUPPRESS_WARNINGS_POP _Pragma("GCC diagnostic pop")
#define PLUTO_SUPPRESS_WARNINGS_TEMPLATE_ID_CDTOR _Pragma("GCC diagnostic ignored \"-Wtemplate-id-cdtor\"")
#define PLUTO_SUPPRESS_WARNINGS_UNUSED_PARAMETER _Pragma("GCC diagnostic ignored \"-Wunused-parameter\"")
#define PLUTO_SUPPRESS_WARNINGS_SIGN_COMPARE _Pragma("GCC diagnostic ignored \"-Wsign-compare\"")
#endif


#if !defined(PLUTO_SUPPRESS_WARNINGS_PUSH)
#define PLUTO_SUPPRESS_WARNINGS_PUSH
#endif
#if !defined(PLUTO_SUPPRESS_WARNINGS_POP)
#define PLUTO_SUPPRESS_WARNINGS_POP
#endif
#if !defined(PLUTO_SUPPRESS_WARNINGS_INTEGER_SIGN_CHANGE)
#define PLUTO_SUPPRESS_WARNINGS_INTEGER_SIGN_CHANGE
#endif
#if !defined(PLUTO_SUPPRESS_WARNINGS_CODE_IS_UNREACHABLE)
#define PLUTO_SUPPRESS_WARNINGS_CODE_IS_UNREACHABLE
#endif
#if !defined(PLUTO_SUPPRESS_WARNINGS_TEMPLATE_ID_CDTOR)
#define PLUTO_SUPPRESS_WARNINGS_TEMPLATE_ID_CDTOR
#endif
#if !defined(PLUTO_SUPPRESS_WARNINGS_UNUSED_PARAMETER)
#define PLUTO_SUPPRESS_WARNINGS_UNUSED_PARAMETER
#endif
#if !defined(PLUTO_SUPPRESS_WARNINGS_SIGN_COMPARE)
#define PLUTO_SUPPRESS_WARNINGS_SIGN_COMPARE
#endif


PLUTO_SUPPRESS_WARNINGS_PUSH
PLUTO_SUPPRESS_WARNINGS_UNUSED_PARAMETER
PLUTO_SUPPRESS_WARNINGS_SIGN_COMPARE
#include "pluto/memory_resource/detail/yakl/YAKL_Gator.h"
PLUTO_SUPPRESS_WARNINGS_POP

// ------------------------------------------------------------------------------------

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>

#include "pluto/alignment.h"
#include "pluto/memory_resource.h"
#include "pluto/stream.h"

namespace pluto {

struct GatorOptions {
    std::size_t initial_size = 0;
    std::size_t grow_size    = 0;
    std::size_t block_size   = 0;
};

class GatorMemoryResource : public memory_resource {
public:
    GatorMemoryResource(): GatorMemoryResource(get_default_resource()) {}

    GatorMemoryResource(memory_resource* upstream): upstream_(upstream) { init(GatorOptions()); }

    GatorMemoryResource(std::unique_ptr<memory_resource>&& upstream):
        owned_upstream_(std::move(upstream)), upstream_(owned_upstream_.get()) {
        init(GatorOptions());
    }

    GatorMemoryResource(const GatorOptions& options): GatorMemoryResource(options, get_default_resource()) {}

    GatorMemoryResource(const GatorOptions& options, memory_resource* upstream): upstream_(upstream) { init(options); }

    GatorMemoryResource(const GatorOptions& options, std::unique_ptr<memory_resource>&& upstream):
        owned_upstream_(std::move(upstream)), upstream_(owned_upstream_.get()) {
        init(options);
    }

    ~GatorMemoryResource() {
        if (gator_) {
            gator_->finalize();
        }
    }

    yakl::Gator& gator() {
        if (not gator_) {
            gator_ = std::make_unique<yakl::Gator>();
            gator_->init(gator_allocate_, gator_deallocate_, gator_zero_, initial_size_, grow_size_, block_size_);
        }
        return *gator_;
    }

protected:
    void init(const GatorOptions& options) {
        gator_allocate_ = [this](std::size_t bytes) {
            scoped_label label("pool_chunk");
            void* ptr              = upstream_->allocate(bytes, alignment_);
            pool_chunk_bytes_[ptr] = bytes;
            return ptr;
        };
        gator_deallocate_ = [this](void* ptr) {
            scoped_label label("pool_chunk");
            return upstream_->deallocate(ptr, pool_chunk_bytes_[ptr], alignment_);
        };
        gator_zero_ = [](void* /*ptr*/, std::size_t /*bytes*/) {};

        // Allocation alignment is guaranteed when alignment_ is Gator's blockSize
        constexpr std::size_t Gb = 1024 * 1024 * 1024;
        block_size_              = (options.block_size ? options.block_size
                                                       : default_alignment());  //Gator default: 16*sizeof(std::size_t) = 64;
        grow_size_               = (options.grow_size ? options.grow_size : 1 * Gb);
        initial_size_            = (options.initial_size ? options.initial_size : 1 * Gb);
        alignment_               = block_size_;
    }

    void* do_allocate(std::size_t bytes, std::size_t alignment) override {
        if (std::remainder(alignment_, alignment) != 0) {
            throw std::runtime_error("alignment " + std::to_string(alignment) +
                                     " not supported. Gator pool block_size (" + std::to_string(block_size_) +
                                     ") must multiple of requested alignment");
        }
        return gator().allocate(bytes);
    }

    void do_deallocate(void* ptr, std::size_t /*bytes*/, std::size_t /*alignment*/) override { gator().free(ptr); }

    bool do_is_equal(const memory_resource& other) const noexcept override {
        if (this == &other) {
            return true;
        }
        return false;
    }

private:
    std::unique_ptr<memory_resource> owned_upstream_;
    memory_resource* upstream_;
    std::size_t alignment_;                          // needs to be bound to allocate/deallocate
    std::map<void*, std::size_t> pool_chunk_bytes_;  // needs to be bound to allocate/deallocate
    std::size_t initial_size_;
    std::size_t grow_size_;
    std::size_t block_size_;
    std::unique_ptr<yakl::Gator> gator_;
    std::function<void*(std::size_t)> gator_allocate_;
    std::function<void(void*)> gator_deallocate_;
    std::function<void(void*, std::size_t)> gator_zero_;
};

}  // namespace pluto
