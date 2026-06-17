/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/**
 * @file relayout_on_host.cc
 * @brief Host implementations for copying Atlas data between blocked and nonblocked layouts.
 *
 * Host relayout works on Atlas views (`atlas::View`/`atlas::ArrayView`) and on mdspan-like
 * views.  The caller is responsible for providing views with the rank, datatype, and shape
 * requirements documented in `atlas/util/relayout.h`.
 */


// Disable bounds checking for internal views in this file,
// as the public API functions already validate shapes and extents before dispatching here.
#define ATLAS_ARRAYVIEW_BOUNDS_CHECKING 0
#ifndef NDEBUG
#define NDEBUG
#endif

#include "atlas/util/relayout.h"

#include <cstdint>
#include <cstdlib>
#include <type_traits>
#include <cstring>

#include "atlas/library/defines.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/runtime/Log.h"
namespace atlas {

namespace {

constexpr const char* relayout_loop_order_env = "ATLAS_RELAYOUT_LOOP_ORDER";
constexpr const char* relayout_nproma_dispatch_env = "ATLAS_RELAYOUT_NPROMA_DISPATCH";
constexpr const char* relayout_blocked_to_blocked_use_memcpy_env = "ATLAS_RELAYOUT_BLOCKED_TO_BLOCKED_USE_MEMCPY";
constexpr const char* relayout_blocked_nonblocked_use_memcpy_env = "ATLAS_RELAYOUT_BLOCKED_NONBLOCKED_USE_MEMCPY";
constexpr const char* relayout_index_operator_env = "ATLAS_RELAYOUT_INDEX_OPERATOR";
constexpr bool blocked_to_blocked_use_memcpy_default = true;
constexpr bool blocked_nonblocked_use_memcpy_default = false;
constexpr bool index_operator_default = false;

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER) || defined(__NVCOMPILER)
#define ATLAS_RELAYOUT_RESTRICT __restrict__
#else
#define ATLAS_RELAYOUT_RESTRICT
#endif

enum class RelayoutLoopOrder {
    nproma_innermost,
    nproma_outermost
};

enum class RelayoutNpromaDispatch {
    static_dispatch,
    runtime,
    runtime_full_blocks
};

RelayoutLoopOrder relayout_loop_order() {
    static RelayoutLoopOrder cached_loop_order = []() {
         const char* loop_order = std::getenv(relayout_loop_order_env);
         if (loop_order && std::strcmp(loop_order, "nproma_outermost") == 0) {
             return RelayoutLoopOrder::nproma_outermost;
         }
         return RelayoutLoopOrder::nproma_innermost;
    }();
    return cached_loop_order;
}

RelayoutNpromaDispatch relayout_nproma_dispatch() {
    static RelayoutNpromaDispatch cached_dispatch = []() {
         const char* dispatch = std::getenv(relayout_nproma_dispatch_env);
         if (dispatch && std::strcmp(dispatch, "runtime") == 0) {
             return RelayoutNpromaDispatch::runtime;
         }
         if (dispatch && std::strcmp(dispatch, "runtime_full_blocks") == 0) {
             return RelayoutNpromaDispatch::runtime_full_blocks;
         }
         return RelayoutNpromaDispatch::static_dispatch;
    }();
    return cached_dispatch;
}

bool relayout_env_flag(const char* name, const bool default_value) {
    const char* value = std::getenv(name);
    if (value == nullptr) {
        return default_value;
    }
    if (std::strcmp(value, "0") == 0 || std::strcmp(value, "false") == 0 || std::strcmp(value, "off") == 0) {
        return false;
    }
    if (std::strcmp(value, "1") == 0 || std::strcmp(value, "true") == 0 || std::strcmp(value, "on") == 0) {
        return true;
    }
    return default_value;
}

bool relayout_blocked_to_blocked_use_memcpy() {
    static bool cached = relayout_env_flag(relayout_blocked_to_blocked_use_memcpy_env, blocked_to_blocked_use_memcpy_default);
    return cached;
}

bool relayout_blocked_nonblocked_use_memcpy() {
    static bool cached = relayout_env_flag(relayout_blocked_nonblocked_use_memcpy_env, blocked_nonblocked_use_memcpy_default);
    return cached;
}

bool relayout_index_operator() {
    static bool cached = relayout_env_flag(relayout_index_operator_env, index_operator_default);
    return cached;
}

// C++17 SFINAE detector to check if a type is a std::mdspan 
template <typename T, typename = void>
struct is_mdspan : std::false_type {};

template <typename T>
struct is_mdspan<T, std::void_t<
    typename T::element_type,
    typename T::data_handle_type,
    typename T::mapping_type,
    decltype(std::declval<T>().data_handle()),
    decltype(std::declval<T>().accessor())
>> : std::true_type {};

// Detects if an instance exposes a .contiguous() member function (like atlas::ArrayView or atlas::LocalView)
template <typename T, typename = std::void_t<>>
struct has_contiguous_member : std::false_type {};

template <typename T>
struct has_contiguous_member<T, std::void_t<
    decltype(std::declval<const T&>().contiguous())
>> : std::is_convertible<decltype(std::declval<const T&>().contiguous()), bool> {};

// Detects if a type has mdspan-style type definitions
template <typename ViewType>
[[nodiscard]] constexpr bool is_layout_always_contiguous() noexcept {
    if constexpr (is_mdspan<ViewType>::value) {
        using mapping_type = typename ViewType::mapping_type;
        return mapping_type::is_always_unique() && mapping_type::is_always_exhaustive();
    } else {
        // Types like atlas::ArrayView resolve contiguity at runtime, 
        // so they do not offer a strict static type guarantee.
        return false;
    }
}

template <typename>
[[maybe_unused]] inline constexpr bool always_false_v = false;

template <typename ViewType,
          typename = std::enable_if_t<is_mdspan<ViewType>::value || 
                                      has_contiguous_member<ViewType>::value>>
[[nodiscard]] constexpr bool is_contiguous(const ViewType& view) noexcept {
    // Branch A: If the type is statically proven contiguous, compile out the rest.
    if constexpr (is_layout_always_contiguous<ViewType>()) {
        return true;
    } 
    // Branch B: Handle standard mdspan fallbacks (dynamic/strided layouts)
    else if constexpr (is_mdspan<ViewType>::value) {
        return view.is_unique() && view.is_exhaustive();
    } 
    // Branch C: Duck-type handle for atlas::ArrayView or identical APIs
    else if constexpr (has_contiguous_member<ViewType>::value) {
        return view.contiguous();
    }
    else {
        static_assert(always_false_v<ViewType>, "Unsupported view type for is_contiguous");
    }
    return false;
}

// C++17 SFINAE detector to check for .data() member function returning a pointer
template <typename T, typename = void>
struct has_member_data : std::false_type {};

template <typename T>
struct is_pointer_type : std::is_pointer<decltype(std::declval<T>().data())> {};

template <typename T>
struct has_member_data<T, std::void_t<decltype(std::declval<T>().data())>> 
    : is_pointer_type<T> {};


// Access raw data pointer for std::mdspan types and types with .data() member function
// Pathway A: For std::mdspan types
template <typename T, typename std::enable_if_t<is_mdspan<std::decay_t<T>>::value, int> = 0>
constexpr auto* get_raw_data(T&& view) {
    auto handle = view.data_handle();
    auto accessor = view.accessor();
    return &accessor.access(handle, 0);
}

// Pathway B: For array::ArrayView types (or any type offering .data())
template <typename T, typename std::enable_if_t<has_member_data<std::decay_t<T>>::value, int> = 0>
constexpr auto* get_raw_data(T&& view) {
    return view.data();
}

// Check that each block slice [dim1..last] is tightly packed, regardless of spacing between blocks.
template <typename Blocked>
[[nodiscard]] bool is_block_contiguous(const Blocked& blocked) {
    const idx_t rank = blocked.rank();
    ATLAS_ASSERT(rank >= 2, "is_block_contiguous requires rank >= 2");
    idx_t expected_stride = 1;
    for (idx_t dim = rank - 1; dim >= 1; --dim) {
        const idx_t stride = static_cast<idx_t>(blocked.stride(dim));
        if (stride != expected_stride) {
            return false;
        }
        expected_stride *= static_cast<idx_t>(blocked.extent(dim));
    }
    return true;
}

template <idx_t nrof_static = 0, typename BlockedValue, typename NonblockedValue>
void copy_nonblocked_to_blocked_rank4_block(BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                            const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                            [[maybe_unused]] const idx_t nrof, // unused if nrof_static != 0
                                            const idx_t nlev,
                                            const idx_t nvar,
                                            const idx_t nproma,
                                            const RelayoutLoopOrder loop_order) {
    const idx_t point_stride = nlev * nvar;
    const idx_t var_block_stride = nlev * nproma;
    const idx_t lev_block_stride = nproma;
    const idx_t lev_nonblocked_stride = nvar;
    if (loop_order == RelayoutLoopOrder::nproma_innermost) {
        for (idx_t jvar = 0; jvar < nvar; ++jvar) {
            BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jvar = raw_blocked_jblk + jvar * var_block_stride;
            const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jvar = raw_nonblocked_jblk + jvar;
            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jvar_jlev = raw_blocked_jblk_jvar + jlev * lev_block_stride;
                const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jvar_jlev = raw_nonblocked_jblk_jvar + jlev * lev_nonblocked_stride;
                if constexpr (nrof_static == 0) {
                    atlas_omp_pragma(omp simd)
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        raw_blocked_jblk_jvar_jlev[jrof] = raw_nonblocked_jblk_jvar_jlev[jrof * point_stride];
                    }
                }
                else {
                    atlas_omp_pragma(omp simd)
                    for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                        raw_blocked_jblk_jvar_jlev[jrof] = raw_nonblocked_jblk_jvar_jlev[jrof * point_stride];
                    }
                }
            }
        }
    }
    else {
        if constexpr (nrof_static == 0) {
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * point_stride;
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof_jvar = raw_blocked_jblk_jrof + jvar * var_block_stride;
                    const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof_jvar = raw_nonblocked_jblk_jrof + jvar;
                    atlas_omp_pragma(omp simd)
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        raw_blocked_jblk_jrof_jvar[jlev * lev_block_stride] = raw_nonblocked_jblk_jrof_jvar[jlev * lev_nonblocked_stride];
                    }
                }
            }
        }
        else {
            for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * point_stride;
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof_jvar = raw_blocked_jblk_jrof + jvar * var_block_stride;
                    const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof_jvar = raw_nonblocked_jblk_jrof + jvar;
                    atlas_omp_pragma(omp simd)
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        raw_blocked_jblk_jrof_jvar[jlev * lev_block_stride] = raw_nonblocked_jblk_jrof_jvar[jlev * lev_nonblocked_stride];
                    }
                }
            }
        }
    }
}

template <idx_t nrof_static = 0, typename BlockedValue, typename NonblockedValue>
void copy_nonblocked_to_blocked_rank3_block(BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                            const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                            [[maybe_unused]] const idx_t nrof, // unused if nrof_static != 0
                                            const idx_t nlev,
                                            [[maybe_unused]] const idx_t nproma, // unused if nrof_static != 0
                                            const RelayoutLoopOrder loop_order) {
    constexpr auto nproma_static = nrof_static; // if nrof_static, then nproma is a compile-time constant equal to nrof_static
    if (loop_order == RelayoutLoopOrder::nproma_innermost) {
        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
            const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jlev = raw_nonblocked_jblk + jlev;
            if constexpr (nrof_static) {
                BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jlev = raw_blocked_jblk + jlev * nproma_static;
                atlas_omp_pragma(omp simd)
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    raw_blocked_jblk_jlev[jrof] = raw_nonblocked_jblk_jlev[jrof * nlev];
                }
            }
            else {
                BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jlev = raw_blocked_jblk + jlev * nproma;
                atlas_omp_pragma(omp simd)
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    raw_blocked_jblk_jlev[jrof] = raw_nonblocked_jblk_jlev[jrof * nlev];
                }
            }
        }
    }
    else {
        if constexpr (nrof_static) {
            for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * nlev;
                atlas_omp_pragma(omp simd)
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    raw_blocked_jblk_jrof[jlev * nproma_static] = raw_nonblocked_jblk_jrof[jlev];
                }
            }
        }
        else {
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * nlev;
                atlas_omp_pragma(omp simd)
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    raw_blocked_jblk_jrof[jlev * nproma] = raw_nonblocked_jblk_jrof[jlev];
                }
            }
        }
    }
}

template <typename BlockedValue, typename NonblockedValue>
void copy_nonblocked_to_blocked_rank4_full_block(BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                                 const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                                 const idx_t nlev,
                                                 const idx_t nvar,
                                                 const idx_t nproma,
                                                 const RelayoutLoopOrder loop_order) {
    // We only get here if nproma_dispatch == RelayoutNpromaDispatch::runtime_full_blocks
    switch (nproma) {
        case 8:    return copy_nonblocked_to_blocked_rank4_block<8>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        case 16:   return copy_nonblocked_to_blocked_rank4_block<16>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        case 32:   return copy_nonblocked_to_blocked_rank4_block<32>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        case 64:   return copy_nonblocked_to_blocked_rank4_block<64>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        case 128:  return copy_nonblocked_to_blocked_rank4_block<128>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        case 256:  return copy_nonblocked_to_blocked_rank4_block<256>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        default:   return copy_nonblocked_to_blocked_rank4_block(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
    }
}

template <typename BlockedValue, typename NonblockedValue>
void copy_nonblocked_to_blocked_rank3_full_block(BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                                 const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                                 const idx_t nlev,
                                                 const idx_t nproma,
                                                 const RelayoutLoopOrder loop_order) {
    // We only get here if nproma_dispatch == RelayoutNpromaDispatch::runtime_full_blocks
    switch (nproma) {
        case 8:    return copy_nonblocked_to_blocked_rank3_block<8>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        case 16:   return copy_nonblocked_to_blocked_rank3_block<16>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        case 32:   return copy_nonblocked_to_blocked_rank3_block<32>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        case 64:   return copy_nonblocked_to_blocked_rank3_block<64>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        case 128:  return copy_nonblocked_to_blocked_rank3_block<128>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        case 256:  return copy_nonblocked_to_blocked_rank3_block<256>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        default:   return copy_nonblocked_to_blocked_rank3_block(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
    }
}

template <idx_t nproma_static = 0, class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_contiguous_mdspan(const Nonblocked nonblocked, Blocked blocked) {
    // Requirements for this implementation are that the blocked view is contiguous in the last rank, and the nonblocked view is contiguous in the first rank.
    ATLAS_ASSERT(is_block_contiguous(blocked));
    ATLAS_ASSERT(is_contiguous(nonblocked));

    const idx_t nproma = blocked.extent(blocked.rank()-1);
    if constexpr (nproma_static != 0) {
        ATLAS_ASSERT(nproma_static == nproma);
    }

    const idx_t np     = nonblocked.extent(0);
    const idx_t nblks  = blocked.extent(0);
    [[maybe_unused]] const RelayoutLoopOrder loop_order = relayout_loop_order();
    [[maybe_unused]] const RelayoutNpromaDispatch nproma_dispatch = relayout_nproma_dispatch();
    [[maybe_unused]] const bool use_memcpy = relayout_blocked_nonblocked_use_memcpy();
    static_assert(nonblocked.rank() == blocked.rank()-1);

    if constexpr(blocked.rank()==4) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
        const idx_t nlev = nonblocked.extent(1);
        const idx_t nvar = nonblocked.extent(2);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            const idx_t jpbegin = jblk * nproma;
            auto* raw_blocked_jblk = &blocked(jblk, 0, 0, 0);
            auto* raw_nonblocked_jblk = &nonblocked(jpbegin, 0, 0);
            if (jblk < nblks-1) {
                if constexpr (nproma_static != 0) {
                    copy_nonblocked_to_blocked_rank4_block<nproma_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma_static, loop_order);
                }
                else {
                    if (nproma_dispatch == RelayoutNpromaDispatch::runtime_full_blocks) {
                        copy_nonblocked_to_blocked_rank4_full_block(
                            raw_blocked_jblk, raw_nonblocked_jblk, nlev, nvar, nproma, loop_order);
                    }
                    else {
                        copy_nonblocked_to_blocked_rank4_block(
                            raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
                    }
                }
            }
            else {
                // Last block may be partial, so compute nrof on the fly to avoid out-of-bounds accesses.
                const idx_t nrof = std::min(np - jpbegin, nproma);
                if constexpr (nproma_static != 0) {
                    copy_nonblocked_to_blocked_rank4_block<nproma_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nvar, nproma_static, loop_order);
                }
                else {
                    copy_nonblocked_to_blocked_rank4_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nvar, nproma, loop_order);
                }
            }
        }
    }
    else if constexpr (blocked.rank()==3) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
        const idx_t nlev = nonblocked.extent(1);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            const idx_t jpbegin = jblk * nproma;
            auto* raw_blocked_jblk = &blocked(jblk, 0, 0);
            auto* raw_nonblocked_jblk = &nonblocked(jpbegin, 0);
            if (jblk < nblks-1) {
                if constexpr (nproma_static != 0) {
                    copy_nonblocked_to_blocked_rank3_block<nproma_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma_static, loop_order);
                }
                else {
                    if (nproma_dispatch == RelayoutNpromaDispatch::runtime_full_blocks) {
                        copy_nonblocked_to_blocked_rank3_full_block(
                            raw_blocked_jblk, raw_nonblocked_jblk, nlev, nproma, loop_order);
                    }
                    else {
                        copy_nonblocked_to_blocked_rank3_block(
                            raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
                    }
                }
            }
            else { // last block may not be complete
                const idx_t nrof = std::min(np - jpbegin, nproma);
                if constexpr (nproma_static != 0) {
                    copy_nonblocked_to_blocked_rank3_block<nproma_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nproma_static, loop_order);
                }
                else {
                    copy_nonblocked_to_blocked_rank3_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nproma, loop_order);
                }
            }
        }
    }
    else if constexpr (blocked.rank()==2) {
        using Value = std::decay_t<decltype(nonblocked(0))>;
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            const idx_t jpbegin = jblk * nproma;
            auto* raw_blocked_jblk = &blocked(jblk, 0);
            auto* raw_nonblocked_jblk = &nonblocked(jpbegin);
            if constexpr (nproma_static != 0) {
                if (use_memcpy) {
                    std::memcpy(
                        raw_blocked_jblk,
                        raw_nonblocked_jblk,
                        static_cast<std::size_t>(nproma_static) * sizeof(Value));
                }
                else {
                    for (idx_t jrof = 0; jrof < nproma_static; ++jrof) {
                        raw_blocked_jblk[jrof] = raw_nonblocked_jblk[jrof];
                    }
                }

            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                if (use_memcpy) {
                    std::memcpy(
                        raw_blocked_jblk,
                        raw_nonblocked_jblk,
                        static_cast<std::size_t>(nrof) * sizeof(Value));
                }
                else {
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        raw_blocked_jblk[jrof] = raw_nonblocked_jblk[jrof];
                    }
                }
            }
        }
    }
    else {
        ATLAS_THROW_EXCEPTION("host_copy_nonblocked_to_blocked_mdspan not implemented");
    }
}

template <idx_t nrof_static = 0, typename BlockedValue, typename NonblockedValue>
void copy_blocked_to_nonblocked_rank4_block(const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                            NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                            [[maybe_unused]] const idx_t nrof, // unused if nrof_static != 0
                                            const idx_t nlev,
                                            const idx_t nvar,
                                            const idx_t nproma,
                                            const RelayoutLoopOrder loop_order) {
    const idx_t point_stride = nlev * nvar;
    const idx_t var_block_stride = nlev * nproma;
    const idx_t lev_block_stride = nproma;
    const idx_t lev_nonblocked_stride = nvar;

    if (loop_order == RelayoutLoopOrder::nproma_innermost) {
        for (idx_t jvar = 0; jvar < nvar; ++jvar) {
            const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jvar = raw_blocked_jblk + jvar * var_block_stride;
            NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jvar = raw_nonblocked_jblk + jvar;
            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jvar_jlev = raw_blocked_jblk_jvar + jlev * lev_block_stride;
                NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jvar_jlev = raw_nonblocked_jblk_jvar + jlev * lev_nonblocked_stride;
                if constexpr (nrof_static == 0) {
                    atlas_omp_pragma(omp simd)
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        raw_nonblocked_jblk_jvar_jlev[jrof * point_stride] = raw_blocked_jblk_jvar_jlev[jrof];
                    }
                }
                else {
                    atlas_omp_pragma(omp simd)
                    for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                        raw_nonblocked_jblk_jvar_jlev[jrof * point_stride] = raw_blocked_jblk_jvar_jlev[jrof];
                    }
                }
            }
        }
    }
    else {
        if constexpr (nrof_static == 0) {
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * point_stride;
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof_jvar = raw_blocked_jblk_jrof + jvar * var_block_stride;
                    NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof_jvar = raw_nonblocked_jblk_jrof + jvar;
                    atlas_omp_pragma(omp simd)
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        raw_nonblocked_jblk_jrof_jvar[jlev * lev_nonblocked_stride] = raw_blocked_jblk_jrof_jvar[jlev * lev_block_stride];
                    }
                }
            }
        }
        else {
            for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * point_stride;
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof_jvar = raw_blocked_jblk_jrof + jvar * var_block_stride;
                    NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof_jvar = raw_nonblocked_jblk_jrof + jvar;
                    atlas_omp_pragma(omp simd)
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        raw_nonblocked_jblk_jrof_jvar[jlev * lev_nonblocked_stride] = raw_blocked_jblk_jrof_jvar[jlev * lev_block_stride];
                    }
                }
            }
        }
    }
}

template <idx_t nrof_static = 0, typename BlockedValue, typename NonblockedValue>
void copy_blocked_to_nonblocked_rank3_block(const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                            NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                            [[maybe_unused]] const idx_t nrof, // unused if nrof_static != 0
                                            const idx_t nlev,
                                            const idx_t nproma,
                                            const RelayoutLoopOrder loop_order) {
    if (loop_order == RelayoutLoopOrder::nproma_innermost) {
        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
            const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jlev = raw_blocked_jblk + jlev * nproma;
            NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jlev = raw_nonblocked_jblk + jlev;
            if constexpr (nrof_static == 0) {
                atlas_omp_pragma(omp simd)
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    raw_nonblocked_jblk_jlev[jrof * nlev] = raw_blocked_jblk_jlev[jrof];
                }
            }
            else {
                atlas_omp_pragma(omp simd)
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    raw_nonblocked_jblk_jlev[jrof * nlev] = raw_blocked_jblk_jlev[jrof];
                }
            }
        }
    }
    else {
        if constexpr (nrof_static == 0) {
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * nlev;
                atlas_omp_pragma(omp simd)
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    raw_nonblocked_jblk_jrof[jlev] = raw_blocked_jblk_jrof[jlev * nproma];
                }
            }
        }
        else {
            for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * nlev;
                atlas_omp_pragma(omp simd)
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    raw_nonblocked_jblk_jrof[jlev] = raw_blocked_jblk_jrof[jlev * nproma];
                }
            }
        }
    }
}

template <typename BlockedValue, typename NonblockedValue>
void copy_blocked_to_nonblocked_rank4_full_block(const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                                 NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                                 const idx_t nlev,
                                                 const idx_t nvar,
                                                 const idx_t nproma,
                                                 const RelayoutLoopOrder loop_order) {
    // We only get here if nproma_dispatch == RelayoutNpromaDispatch::runtime_full_blocks
    switch (nproma) {
        case 8:    return copy_blocked_to_nonblocked_rank4_block<8>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        case 16:   return copy_blocked_to_nonblocked_rank4_block<16>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        case 32:   return copy_blocked_to_nonblocked_rank4_block<32>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        case 64:   return copy_blocked_to_nonblocked_rank4_block<64>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        case 128:  return copy_blocked_to_nonblocked_rank4_block<128>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        case 256:  return copy_blocked_to_nonblocked_rank4_block<256>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
        default:   return copy_blocked_to_nonblocked_rank4_block(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
    }
}

template <typename BlockedValue, typename NonblockedValue>
void copy_blocked_to_nonblocked_rank3_full_block(const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                                 NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                                 const idx_t nlev,
                                                 const idx_t nproma,
                                                 const RelayoutLoopOrder loop_order) {
    // We only get here if nproma_dispatch == RelayoutNpromaDispatch::runtime_full_blocks
    switch (nproma) {
        case 8:    return copy_blocked_to_nonblocked_rank3_block<8>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        case 16:   return copy_blocked_to_nonblocked_rank3_block<16>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        case 32:   return copy_blocked_to_nonblocked_rank3_block<32>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        case 64:   return copy_blocked_to_nonblocked_rank3_block<64>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        case 128:  return copy_blocked_to_nonblocked_rank3_block<128>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        case 256:  return copy_blocked_to_nonblocked_rank3_block<256>(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
        default:   return copy_blocked_to_nonblocked_rank3_block(raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
    }
}

template <idx_t nproma_static = 0, class Blocked, class Nonblocked>
void host_copy_blocked_to_nonblocked_contiguous_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    // Requirements: blocked and nonblocked must be contiguous, and if nproma_static is specified, it must match the last dimension of blocked.
    ATLAS_ASSERT(is_block_contiguous(blocked));
    ATLAS_ASSERT(is_contiguous(nonblocked));

    const idx_t nproma = blocked.extent(blocked.rank()-1);
    if constexpr (nproma_static != 0) {
        ATLAS_ASSERT(nproma_static == nproma);
    }

    const idx_t np     = nonblocked.extent(0);
    const idx_t nblks  = blocked.extent(0);
    [[maybe_unused]] const RelayoutLoopOrder loop_order = relayout_loop_order();
    [[maybe_unused]] const RelayoutNpromaDispatch nproma_dispatch = relayout_nproma_dispatch();
    [[maybe_unused]] const bool use_memcpy = relayout_blocked_nonblocked_use_memcpy();
    static_assert(nonblocked.rank() == blocked.rank()-1);

    if constexpr(blocked.rank()==4) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
        const idx_t nlev = nonblocked.extent(1);
        const idx_t nvar = nonblocked.extent(2);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            const idx_t jpbegin = jblk * nproma;
            const auto* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk = &blocked(jblk, 0, 0, 0);
            auto* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk = &nonblocked(jpbegin, 0, 0);
            if (jblk < nblks-1) {
                if constexpr (nproma_static != 0) {
                    copy_blocked_to_nonblocked_rank4_block<nproma_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma_static, loop_order);
                }
                else if (nproma_dispatch == RelayoutNpromaDispatch::runtime_full_blocks) {
                    copy_blocked_to_nonblocked_rank4_full_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nlev, nvar, nproma, loop_order);
                }
                else {
                    copy_blocked_to_nonblocked_rank4_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
                }
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                if constexpr (nproma_static != 0) {
                    copy_blocked_to_nonblocked_rank4_block<nproma_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nvar, nproma_static, loop_order);
                }
                else {
                    copy_blocked_to_nonblocked_rank4_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nvar, nproma, loop_order);
                }
            }
        }
    }
    else if constexpr (blocked.rank()==3) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
        const idx_t nlev = nonblocked.extent(1);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            const idx_t jpbegin = jblk * nproma;
            const auto* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk = &blocked(jblk, 0, 0);
            auto* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk = &nonblocked(jpbegin, 0);
            if (jblk < nblks-1) {
                if constexpr (nproma_static != 0) {
                    copy_blocked_to_nonblocked_rank3_block<nproma_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma_static, loop_order);
                }
                else if (nproma_dispatch == RelayoutNpromaDispatch::runtime_full_blocks) {
                    copy_blocked_to_nonblocked_rank3_full_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nlev, nproma, loop_order);
                }
                else {
                    copy_blocked_to_nonblocked_rank3_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
                }
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                if constexpr (nproma_static != 0) {
                    copy_blocked_to_nonblocked_rank3_block<nproma_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nproma_static, loop_order);
                }
                else {
                    copy_blocked_to_nonblocked_rank3_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nproma, loop_order);
                }
            }
        }
    }
    else if constexpr (blocked.rank()==2) {
        using Value = std::remove_cv_t<std::remove_reference_t<decltype(blocked(0, 0))>>;
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            const idx_t jpbegin = jblk * nproma;
            const idx_t nrof = std::min(np - jpbegin, nproma);
            const auto* raw_blocked_jblk = &blocked(jblk, 0);
            auto* raw_nonblocked_jblk = &nonblocked(jpbegin);
            if (use_memcpy) {
                std::memcpy(
                    raw_nonblocked_jblk,
                    raw_blocked_jblk,
                    static_cast<std::size_t>(nrof) * sizeof(Value));
            }
            else {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    raw_nonblocked_jblk[jrof] = raw_blocked_jblk[jrof];
                }
            }
        }
    }
    else {
        ATLAS_THROW_EXCEPTION("host_copy_blocked_to_nonblocked_mdspan not implemented");
    }
}

template <idx_t nproma, class Blocked, class Nonblocked>
void host_copy_blocked_to_nonblocked_contiguous_mdspan_nproma(const Blocked blocked, Nonblocked nonblocked) {
    host_copy_blocked_to_nonblocked_contiguous_mdspan<nproma>(blocked, nonblocked);
}

}  // namespace


/**
 * @brief Copy a nonblocked host view into a blocked host view.
 *
 * @param nonblocked Source Atlas view (`atlas::View`/`atlas::ArrayView`) or mdspan-like view
 *        with rank one less than `blocked`.  For rank-3 views, the dimension order is
 *        `[npoint, nlev, nvar]`.
 * @param blocked Target Atlas view (`atlas::View`/`atlas::ArrayView`) or mdspan-like view with
 *        rank 2, 3, or 4.  For rank-4 views, the dimension order is
 *        `[nblk, nvar, nlev, nproma]`.
 *
 * @pre `nonblocked.rank() == blocked.rank() - 1`.
 * @pre Shared dimensions must match; for rank 4, `nonblocked.extent(1) == blocked.extent(2)`
 *      and `nonblocked.extent(2) == blocked.extent(1)`.
 */
template <idx_t nproma, class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_contiguous_mdspan_nproma(const Nonblocked nonblocked, Blocked blocked) {
    host_copy_nonblocked_to_blocked_contiguous_mdspan<nproma>(nonblocked, blocked);
}

template <class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked) {
    const idx_t nproma = blocked.extent(blocked.rank()-1);
    const idx_t np     = nonblocked.extent(0);
    const idx_t nblks  = blocked.extent(0);
    static_assert(nonblocked.rank() == blocked.rank()-1);

    if (!relayout_index_operator() && is_contiguous(nonblocked) && is_block_contiguous(blocked)) {
        if (relayout_nproma_dispatch() != RelayoutNpromaDispatch::static_dispatch) {
            return host_copy_nonblocked_to_blocked_contiguous_mdspan(nonblocked, blocked);
        }

        // Optimized paths possible with static nproma dispatch.
        switch (nproma) {
            case 8:   return host_copy_nonblocked_to_blocked_contiguous_mdspan_nproma<8 >(nonblocked, blocked);
            case 16:  return host_copy_nonblocked_to_blocked_contiguous_mdspan_nproma<16>(nonblocked, blocked);
            case 32:  return host_copy_nonblocked_to_blocked_contiguous_mdspan_nproma<32>(nonblocked, blocked);
            case 64:  return host_copy_nonblocked_to_blocked_contiguous_mdspan_nproma<64>(nonblocked, blocked);
            case 128: return host_copy_nonblocked_to_blocked_contiguous_mdspan_nproma<128>(nonblocked, blocked);
            case 256: return host_copy_nonblocked_to_blocked_contiguous_mdspan_nproma<256>(nonblocked, blocked);
            default:  return host_copy_nonblocked_to_blocked_contiguous_mdspan(nonblocked, blocked);
        }
    }

    // Fall back to generic implementation if the views are not contiguous or block-contiguous.
    [[maybe_unused]] const RelayoutLoopOrder loop_order = relayout_loop_order();

    if constexpr(blocked.rank()==4) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);

        if (loop_order == RelayoutLoopOrder::nproma_outermost) {
            atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
                const idx_t jpbegin = jblk * nproma;
                const idx_t nrof = std::min(np - jpbegin, nproma);
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        for (idx_t jrof = 0, jp = jpbegin; jrof < nrof; ++jrof, ++jp) {
                            blocked(jblk, jvar, jlev, jrof) = nonblocked(jp, jlev, jvar);
                        }
                    }
                }
            }
        } else {
            atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
                const idx_t jpbegin = jblk * nproma;
                const idx_t nrof = std::min(np - jpbegin, nproma);
                for (idx_t jrof = 0, jp = jpbegin; jrof < nrof; ++jrof, ++jp) {
                    for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            blocked(jblk, jvar, jlev, jrof) = nonblocked(jp, jlev, jvar);
                        }
                    }
                }
            }
        }
    }
    else if constexpr (blocked.rank()==3) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
        idx_t nlev = nonblocked.extent(1);

        if (loop_order == RelayoutLoopOrder::nproma_outermost) {
            atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
                const idx_t jpbegin = jblk * nproma;
                const idx_t nrof = std::min(np - jpbegin, nproma);
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    for (idx_t jrof = 0, jp = jpbegin; jrof < nrof; ++jrof, ++jp) {
                        blocked(jblk, jlev, jrof) = nonblocked(jp, jlev);
                    }
                }
            }
        } else {
            atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
                const idx_t jpbegin = jblk * nproma;
                const idx_t nrof = std::min(np - jpbegin, nproma);
                for (idx_t jrof = 0, jp = jpbegin; jrof < nrof; ++jrof, ++jp) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        blocked(jblk, jlev, jrof) = nonblocked(jp, jlev);
                    }
                }
            }
        }
    }
    else if constexpr (blocked.rank()==2) {
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            const idx_t jpbegin = jblk * nproma;
            const idx_t nrof = std::min(np - jpbegin, nproma);
            for (idx_t jrof = 0, jp = jpbegin; jrof < nrof; ++jrof, ++jp) {
                blocked(jblk, jrof) = nonblocked(jp);
            }
        }
    }
    else {
        ATLAS_THROW_EXCEPTION("host_copy_blocked_to_nonblocked_mdspan not implemented");
    }
}


template <class Blocked, class Nonblocked>
/**
 * @brief Copy a blocked host view into a nonblocked host view.
 *
 * @param blocked Source Atlas view (`atlas::View`/`atlas::ArrayView`) or mdspan-like view with
 *        rank 2, 3, or 4.  For rank-4 views, the dimension order is
 *        `[nblk, nvar, nlev, nproma]`.
 * @param nonblocked Target Atlas view (`atlas::View`/`atlas::ArrayView`) or mdspan-like view
 *        with rank one less than `blocked`.  For rank-3 views, the dimension order is
 *        `[npoint, nlev, nvar]`.
 *
 * @pre `nonblocked.rank() == blocked.rank() - 1`.
 * @pre Shared dimensions must match; for rank 4, `nonblocked.extent(1) == blocked.extent(2)`
 *      and `nonblocked.extent(2) == blocked.extent(1)`.
 */
void host_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    auto np     = nonblocked.extent(0);
    auto nblks  = blocked.extent(0);
    auto nproma = blocked.extent(blocked.rank()-1);
    static_assert(nonblocked.rank() == blocked.rank()-1);

    if (!relayout_index_operator() && is_contiguous(nonblocked) && is_block_contiguous(blocked)) {
        if (relayout_nproma_dispatch() != RelayoutNpromaDispatch::static_dispatch) {
            return host_copy_blocked_to_nonblocked_contiguous_mdspan(blocked, nonblocked);
        }

        // Optimized paths possible with static nproma dispatch.
        switch (nproma) {
            case 8:   return host_copy_blocked_to_nonblocked_contiguous_mdspan_nproma<8  >(blocked, nonblocked);
            case 16:  return host_copy_blocked_to_nonblocked_contiguous_mdspan_nproma<16 >(blocked, nonblocked);
            case 32:  return host_copy_blocked_to_nonblocked_contiguous_mdspan_nproma<32 >(blocked, nonblocked);
            case 64:  return host_copy_blocked_to_nonblocked_contiguous_mdspan_nproma<64 >(blocked, nonblocked);
            case 128: return host_copy_blocked_to_nonblocked_contiguous_mdspan_nproma<128>(blocked, nonblocked);
            case 256: return host_copy_blocked_to_nonblocked_contiguous_mdspan_nproma<256>(blocked, nonblocked);
            default:  return host_copy_blocked_to_nonblocked_contiguous_mdspan(blocked, nonblocked);
        }
    }

    // Fall back to generic implementation if the views are not contiguous or block-contiguous.
    [[maybe_unused]] const RelayoutLoopOrder loop_order = relayout_loop_order();

    if constexpr(blocked.rank()==4) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);

        if (loop_order == RelayoutLoopOrder::nproma_outermost) {
            atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
                const idx_t jpbegin = jblk * nproma;
                const idx_t nrof = std::min(np - jpbegin, nproma);
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        for (idx_t jrof = 0, jp = jpbegin; jrof < nrof; ++jrof, ++jp) {
                            nonblocked(jp, jlev, jvar) = blocked(jblk, jvar, jlev, jrof);
                        }
                    }
                }
            }
        } else {
            atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
                const idx_t jpbegin = jblk * nproma;
                const idx_t nrof = std::min(np - jpbegin, nproma);
                for (idx_t jrof = 0, jp = jpbegin; jrof < nrof; ++jrof, ++jp) {
                    for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            nonblocked(jp, jlev, jvar) = blocked(jblk, jvar, jlev, jrof);
                        }
                    }
                }
            }
        }
    }
    else if constexpr (blocked.rank()==3) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
        idx_t nlev = nonblocked.extent(1);

        if (loop_order == RelayoutLoopOrder::nproma_outermost) {
            atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
                const idx_t jpbegin = jblk * nproma;
                const idx_t nrof = std::min(np - jpbegin, nproma);
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    for (idx_t jrof = 0, jp = jpbegin; jrof < nrof; ++jrof, ++jp) {
                        nonblocked(jp, jlev) = blocked(jblk, jlev, jrof);
                    }
                }
            }
        } else {
            atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
                const idx_t jpbegin = jblk * nproma;
                const idx_t nrof = std::min(np - jpbegin, nproma);
                for (idx_t jrof = 0, jp = jpbegin; jrof < nrof; ++jrof, ++jp) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        nonblocked(jp, jlev) = blocked(jblk, jlev, jrof);
                    }
                }
            }
        }
    }
    else if constexpr (blocked.rank()==2) {
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            const idx_t jpbegin = jblk * nproma;
            const idx_t nrof = std::min(np - jpbegin, nproma);
            for (idx_t jrof = 0, jp = jpbegin; jrof < nrof; ++jrof, ++jp) {
                nonblocked(jp) = blocked(jblk, jrof);
            }
        }
    }
    else {
        ATLAS_THROW_EXCEPTION("host_copy_blocked_to_nonblocked_mdspan not implemented");
    }
}

template <class BlockedIn, class BlockedOut>
/**
 * @brief Copy between two blocked host views, allowing different `nproma` values.
 *
 * @param blocked_in Source Atlas view (`atlas::View`/`atlas::ArrayView`) or mdspan-like blocked
 *        view with rank 2, 3, or 4.
 * @param blocked_out Target Atlas view (`atlas::View`/`atlas::ArrayView`) or mdspan-like blocked
 *        view with the same rank and value type as `blocked_in`.
 *
 * @pre `blocked_in.rank() == blocked_out.rank()`.
 * @pre The value types of the two views must match.
 * @pre Non-horizontal dimensions must match; for rank 4 this means matching `nvar` and `nlev`.
 */
void host_copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out) {
    static_assert(std::is_same_v<std::decay_t<typename BlockedIn::value_type>, std::decay_t<typename BlockedOut::value_type>>, "Data types of input and output views must match for blocked-to-blocked copy");
    using Value = std::decay_t<typename BlockedOut::value_type>;
    const idx_t nblks_in  = blocked_in.extent(0);
    const idx_t nproma_in = blocked_in.extent(blocked_in.rank()-1);
    const idx_t nblks_out  = blocked_out.extent(0);
    const idx_t nproma_out = blocked_out.extent(blocked_out.rank()-1);
    const idx_t total_points_in  = nblks_in * nproma_in;
    const idx_t total_points_out = nblks_out * nproma_out;
    const idx_t total_points     = std::min(total_points_in, total_points_out);
    static_assert(blocked_in.rank() == blocked_out.rank());

    if (total_points == 0) {
        return;
    }

    const bool use_memcpy = relayout_blocked_to_blocked_use_memcpy();

    if (use_memcpy) {
        if (nproma_in == nproma_out && blocked_in.size() == blocked_out.size() && is_contiguous(blocked_in) && is_contiguous(blocked_out)) {
            const Value* raw_in = get_raw_data(blocked_in);
            Value* raw_out = get_raw_data(blocked_out);
            std::memcpy(raw_out, raw_in, blocked_out.size() * sizeof(Value));
            return;
        }
    }

    // At the moment we implement only contiguous block copies for optimizations
    ATLAS_ASSERT(is_block_contiguous(blocked_in));
    ATLAS_ASSERT(is_block_contiguous(blocked_out));

    if constexpr (blocked_in.rank()==4) {
        ATLAS_ASSERT(blocked_in.extent(1) == blocked_out.extent(1));
        ATLAS_ASSERT(blocked_in.extent(2) == blocked_out.extent(2));
        const idx_t nlev = blocked_in.extent(1);
        const idx_t nvar = blocked_in.extent(2);
        const idx_t out_lev_stride = nvar * nproma_out;
        const idx_t out_var_stride = nproma_out;
        const idx_t in_lev_stride = nvar * nproma_in;
        const idx_t in_var_stride = nproma_in;

        atlas_omp_parallel_for(idx_t jblk_out = 0; jblk_out < nblks_out; ++jblk_out) {
            const idx_t jpbegin = jblk_out * nproma_out;
            if (jpbegin >= total_points) {
                continue;
            }

            Value* ATLAS_RELAYOUT_RESTRICT raw_blocked_out_jblk = &blocked_out(jblk_out, 0, 0, 0);

            const idx_t jpend = std::min(total_points, jpbegin + nproma_out);
            idx_t jp = jpbegin;
            while (jp < jpend) {
                const idx_t jblk_in  = jp / nproma_in;
                const idx_t jrof_in  = jp - jblk_in * nproma_in;
                const idx_t jrof_out = jp - jblk_out * nproma_out;
                const idx_t chunk = std::min(jpend - jp, nproma_in - jrof_in);

                const Value* raw_blocked_in_jblk = &blocked_in(jblk_in, 0, 0, 0);

                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                        idx_t index_out_base = jlev * out_lev_stride + jvar * out_var_stride + jrof_out;
                        idx_t index_in_base = jlev * in_lev_stride + jvar * in_var_stride + jrof_in;
                        if (use_memcpy) {
                            std::memcpy(
                                raw_blocked_out_jblk + index_out_base,
                                raw_blocked_in_jblk  + index_in_base,
                                static_cast<std::size_t>(chunk) * sizeof(Value));
                        }
                        else {
                            for (idx_t j = 0; j < chunk; ++j) {
                                idx_t index_out = index_out_base + j;
                                idx_t index_in  = index_in_base  + j;
                                raw_blocked_out_jblk[index_out] = raw_blocked_in_jblk[index_in];
                            }
                        }
                    }
                }
                jp += chunk;
            }
        }
    }
    else if constexpr (blocked_in.rank()==3) {
        ATLAS_ASSERT(blocked_in.extent(1) == blocked_out.extent(1));
        const idx_t nlev = blocked_in.extent(1);
        const idx_t out_lev_stride = nproma_out;
        const idx_t in_lev_stride = nproma_in;

        atlas_omp_parallel_for(idx_t jblk_out = 0; jblk_out < nblks_out; ++jblk_out) {
            const idx_t jpbegin = jblk_out * nproma_out;
            if (jpbegin >= total_points) {
                continue;
            }

            Value* ATLAS_RELAYOUT_RESTRICT raw_blocked_out_jblk = &blocked_out(jblk_out, 0, 0);

            const idx_t jpend = std::min(total_points, jpbegin + nproma_out);
            idx_t jp = jpbegin;
            while (jp < jpend) {
                const idx_t jblk_in  = jp / nproma_in;
                const idx_t jrof_in  = jp - jblk_in * nproma_in;
                const idx_t jrof_out = jp - jblk_out * nproma_out;
                const idx_t chunk = std::min(jpend - jp, nproma_in - jrof_in);

                const Value* ATLAS_RELAYOUT_RESTRICT raw_blocked_in_jblk = &blocked_in(jblk_in, 0, 0);

                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    const idx_t index_out_base = jlev * out_lev_stride + jrof_out;
                    const idx_t index_in_base = jlev * in_lev_stride + jrof_in;
                    if (use_memcpy) {
                        std::memcpy(
                            raw_blocked_out_jblk + index_out_base,
                            raw_blocked_in_jblk + index_in_base,
                            static_cast<std::size_t>(chunk) * sizeof(Value));
                    }
                    else {
                        for (idx_t j = 0; j < chunk; ++j) {
                            raw_blocked_out_jblk[index_out_base + j] = raw_blocked_in_jblk[index_in_base + j];
                        }
                    }
                }
                jp += chunk;
            }
        }
    }
    else if constexpr (blocked_in.rank()==2) {
        atlas_omp_parallel_for(idx_t jblk_out = 0; jblk_out < nblks_out; ++jblk_out) {
            const idx_t jpbegin = jblk_out * nproma_out;
            if (jpbegin >= total_points) {
                continue;
            }

            Value* ATLAS_RELAYOUT_RESTRICT raw_blocked_out_jblk = &blocked_out(jblk_out, 0);

            const idx_t jpend = std::min(total_points, jpbegin + nproma_out);
            idx_t jp = jpbegin;
            while (jp < jpend) {
                const idx_t jblk_in  = jp / nproma_in;
                const idx_t jrof_in  = jp - jblk_in * nproma_in;
                const idx_t jrof_out = jp - jblk_out * nproma_out;
                const idx_t chunk = std::min(jpend - jp, nproma_in - jrof_in);

                const Value* ATLAS_RELAYOUT_RESTRICT raw_blocked_in_jblk = &blocked_in(jblk_in, 0);

                if (use_memcpy) {
                    std::memcpy(
                        raw_blocked_out_jblk + jrof_out,
                        raw_blocked_in_jblk + jrof_in,
                        static_cast<std::size_t>(chunk) * sizeof(Value));
                }
                else {
                    for (idx_t j = 0; j < chunk; ++j) {
                        raw_blocked_out_jblk[jrof_out + j] = raw_blocked_in_jblk[jrof_in + j];
                    }
                }
                jp += chunk;
            }
        }
    }
    else {
        ATLAS_THROW_EXCEPTION("transposition not implemented for rank " << blocked_in.rank());
    }
}

#define EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(TYPE, BLOCKED_RANK) \
    template void host_copy_blocked_to_nonblocked_mdspan<array::ArrayView<const TYPE,BLOCKED_RANK>,array::ArrayView<TYPE,BLOCKED_RANK-1>>(array::ArrayView<const TYPE,BLOCKED_RANK>, array::ArrayView<TYPE,BLOCKED_RANK-1>); \
    template void host_copy_blocked_to_nonblocked_mdspan<mdspan<const TYPE,dims<BLOCKED_RANK>,layout_stride>,mdspan<TYPE,dims<BLOCKED_RANK-1>,layout_stride>>(mdspan<const TYPE,dims<BLOCKED_RANK>,layout_stride>, mdspan<TYPE,dims<BLOCKED_RANK-1>,layout_stride>); \
\
    template void host_copy_nonblocked_to_blocked_mdspan<array::ArrayView<const TYPE,BLOCKED_RANK-1>,array::ArrayView<TYPE,BLOCKED_RANK>>(array::ArrayView<const TYPE,BLOCKED_RANK-1>, array::ArrayView<TYPE,BLOCKED_RANK>); \
    template void host_copy_nonblocked_to_blocked_mdspan<mdspan<const TYPE,dims<BLOCKED_RANK-1>,layout_stride>,mdspan<TYPE,dims<BLOCKED_RANK>,layout_stride>>(mdspan<const TYPE,dims<BLOCKED_RANK-1>,layout_stride>, mdspan<TYPE,dims<BLOCKED_RANK>,layout_stride>); \
\
    template void host_copy_blocked_to_blocked_mdspan<array::ArrayView<const TYPE,BLOCKED_RANK>,array::ArrayView<TYPE,BLOCKED_RANK>>(array::ArrayView<const TYPE,BLOCKED_RANK>, array::ArrayView<TYPE,BLOCKED_RANK>); \
    template void host_copy_blocked_to_blocked_mdspan<mdspan<const TYPE,dims<BLOCKED_RANK>,layout_stride>,mdspan<TYPE,dims<BLOCKED_RANK>,layout_stride>>(mdspan<const TYPE,dims<BLOCKED_RANK>,layout_stride>, mdspan<TYPE,dims<BLOCKED_RANK>,layout_stride>); \

#define EXPLICIT_TEMPLATE_INSTANTIATION(RANK)                \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(double, RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(float , RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(int   , RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(long  , RANK)

EXPLICIT_TEMPLATE_INSTANTIATION(2)
EXPLICIT_TEMPLATE_INSTANTIATION(3)
EXPLICIT_TEMPLATE_INSTANTIATION(4)

#undef EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK
#undef EXPLICIT_TEMPLATE_INSTANTIATION
#undef ATLAS_RELAYOUT_RESTRICT

}  // namespace atlas
