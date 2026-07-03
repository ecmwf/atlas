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
 * @file relayout_on_host.tcc
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

#define DISABLE_RAW_POINTERS 0
#define ATLAS_RELAYOUT_SIMD atlas_omp_pragma(omp simd)
// #define ATLAS_RELAYOUT_SIMD

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
#include "atlas/mdspan.h"
#include "atlas/runtime/Log.h"


#ifdef atlas_omp_parallel_for
#undef atlas_omp_parallel_for
#endif
#define atlas_omp_parallel_for for
namespace atlas {

namespace {

constexpr std::size_t alignment = 64;

constexpr std::size_t greatest_common_divisor(std::size_t lhs, std::size_t rhs) {
    while (rhs != 0) {
        const std::size_t remainder = lhs % rhs;
        lhs = rhs;
        rhs = remainder;
    }
    return lhs;
}

template <size_t nproma, typename Value>
inline constexpr std::size_t blocked_subspan_alignment_v =
    greatest_common_divisor(alignment, static_cast<std::size_t>(nproma) * sizeof(Value));

enum class BlockAlignment { aligned, unaligned };

constexpr const char* relayout_loop_order_env = "ATLAS_RELAYOUT_LOOP_ORDER";
constexpr const char* relayout_nproma_dispatch_env = "ATLAS_RELAYOUT_NPROMA_DISPATCH";
constexpr const char* relayout_blocked_to_blocked_use_memcpy_env = "ATLAS_RELAYOUT_BLOCKED_TO_BLOCKED_USE_MEMCPY";
constexpr const char* relayout_blocked_nonblocked_use_memcpy_env = "ATLAS_RELAYOUT_BLOCKED_NONBLOCKED_USE_MEMCPY";
constexpr const char* relayout_implementation_env = "ATLAS_RELAYOUT_IMPLEMENTATION";
constexpr bool blocked_to_blocked_use_memcpy_default = true;
constexpr bool blocked_nonblocked_use_memcpy_default = false;
constexpr bool index_operator_default = false;
constexpr const char* implementation_default = "raw_pointers";

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
    dynamic,
};

enum class RelayoutImplementation {
    raw_pointers,
    mdspan
};

[[maybe_unused]] RelayoutLoopOrder relayout_loop_order() {
    static RelayoutLoopOrder cached_loop_order = []() {
         const char* loop_order = std::getenv(relayout_loop_order_env);
         if (loop_order && std::strcmp(loop_order, "nproma_outermost") == 0) {
             return RelayoutLoopOrder::nproma_outermost;
         }
         return RelayoutLoopOrder::nproma_innermost;
    }();
    return cached_loop_order;
}

[[maybe_unused]] RelayoutNpromaDispatch relayout_nproma_dispatch() {
    static RelayoutNpromaDispatch cached_dispatch = []() {
         const char* dispatch = std::getenv(relayout_nproma_dispatch_env);
         if (dispatch) {
            if (std::strcmp(dispatch, "static") == 0) {
                return RelayoutNpromaDispatch::static_dispatch;
            }
            else if (std::strcmp(dispatch, "dynamic") == 0) {
                return RelayoutNpromaDispatch::dynamic;
            }
            else {
                ATLAS_THROW_EXCEPTION("Unknown value for ATLAS_RELAYOUT_NPROMA_DISPATCH: " << dispatch
                                    << ".  Accepted values are 'static' or 'dynamic'.");
            }
         }
         return RelayoutNpromaDispatch::static_dispatch;
    }();
    return cached_dispatch;
}

[[maybe_unused]] RelayoutImplementation relayout_implementation() {
    static RelayoutImplementation cached = []() {
         const char* impl = std::getenv(relayout_implementation_env);
         if (impl && std::strcmp(impl, "mdspan") == 0) {
             return RelayoutImplementation::mdspan;
         }
         return RelayoutImplementation::raw_pointers;
    }();
    return cached;
}

[[maybe_unused]] bool relayout_env_flag(const char* name, const bool default_value) {
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

[[maybe_unused]] bool relayout_blocked_to_blocked_use_memcpy() {
    static bool cached = relayout_env_flag(relayout_blocked_to_blocked_use_memcpy_env, blocked_to_blocked_use_memcpy_default);
    return cached;
}

[[maybe_unused]] bool relayout_blocked_nonblocked_use_memcpy() {
    static bool cached = relayout_env_flag(relayout_blocked_nonblocked_use_memcpy_env, blocked_nonblocked_use_memcpy_default);
    return cached;
}

enum class ViewType { Blocked, NonBlocked };

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

// Helper to detect if a view has stride-1 in the last dimension
// This allows it to be treated as layout_right for better performance
template <typename ViewType>
[[nodiscard]] bool has_stride_one_last_dimension(const ViewType& view) {
    const idx_t rank = view.rank();
    if (rank == 0) return true;  // scalar edge case
    return view.stride(rank - 1) == 1;
}

// Check if a view has layout_right-like properties: contiguous and stride-1 in last dimension
// Works for both ArrayView and mdspan types
template <typename ViewType>
[[nodiscard]] bool has_layout_right(const ViewType& view) {
    return is_contiguous(view) && has_stride_one_last_dimension(view);
}

#if DISABLE_RAW_POINTERS == 0
template <idx_t nrof_static = 0, typename BlockedValue, typename NonblockedValue>
void copy_nonblocked_to_blocked_contiguous_rank4_block(BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
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
                    ATLAS_RELAYOUT_SIMD
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        raw_blocked_jblk_jvar_jlev[jrof] = raw_nonblocked_jblk_jvar_jlev[jrof * point_stride];
                    }
                }
                else {
                    ATLAS_RELAYOUT_SIMD
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
                    ATLAS_RELAYOUT_SIMD
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
                    ATLAS_RELAYOUT_SIMD
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        raw_blocked_jblk_jrof_jvar[jlev * lev_block_stride] = raw_nonblocked_jblk_jrof_jvar[jlev * lev_nonblocked_stride];
                    }
                }
            }
        }
    }
}

template <idx_t nrof_static = 0, typename BlockedValue, typename NonblockedValue>
void copy_nonblocked_to_blocked_contiguous_rank3_block(BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
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
                ATLAS_RELAYOUT_SIMD
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    raw_blocked_jblk_jlev[jrof] = raw_nonblocked_jblk_jlev[jrof * nlev];
                }
            }
            else {
                BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jlev = raw_blocked_jblk + jlev * nproma;
                ATLAS_RELAYOUT_SIMD
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
                ATLAS_RELAYOUT_SIMD
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    raw_blocked_jblk_jrof[jlev * nproma_static] = raw_nonblocked_jblk_jrof[jlev];
                }
            }
        }
        else {
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                const NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * nlev;
                ATLAS_RELAYOUT_SIMD
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    raw_blocked_jblk_jrof[jlev * nproma] = raw_nonblocked_jblk_jrof[jlev];
                }
            }
        }
    }
}
#endif

#if DISABLE_RAW_POINTERS == 0
template <size_t nproma_extent, class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma(const Nonblocked nonblocked, Blocked blocked) {
    // Requirements for this implementation are that the blocked view is contiguous in the last rank, and the nonblocked view is contiguous in the first rank.
    ATLAS_ASSERT(is_block_contiguous(blocked));
    ATLAS_ASSERT(has_layout_right(nonblocked));

    const idx_t nproma = blocked.extent(blocked.rank()-1);
    if constexpr (nproma_extent != dynamic_extent) {
        ATLAS_ASSERT(nproma_extent == nproma);
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
                const idx_t nrof = nproma;
                if constexpr (nproma_extent != dynamic_extent) {
                    constexpr idx_t nrof_static = nproma_extent;
                    copy_nonblocked_to_blocked_contiguous_rank4_block<nrof_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nvar, nproma, loop_order);
                }
                else {
                    copy_nonblocked_to_blocked_contiguous_rank4_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nvar, nproma, loop_order);
                }
            }
            else {
                // Last block may be partial, so compute nrof on the fly to avoid out-of-bounds accesses.
                const idx_t nrof = std::min(np - jpbegin, nproma);
                if constexpr (nproma_extent != dynamic_extent) {
                    constexpr idx_t nrof_static = nproma_extent;
                    copy_nonblocked_to_blocked_contiguous_rank4_block<nrof_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nvar, nproma, loop_order);
                }
                else {
                    copy_nonblocked_to_blocked_contiguous_rank4_block(
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
                const idx_t nrof = nproma;
                if constexpr (nproma_extent != dynamic_extent) {
                    constexpr idx_t nrof_static = nproma_extent;
                    copy_nonblocked_to_blocked_contiguous_rank3_block<nrof_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nproma, loop_order);
                }
                else {
                    copy_nonblocked_to_blocked_contiguous_rank3_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nproma, loop_order);
                }
            }
            else { // last block may not be complete
                const idx_t nrof = std::min(np - jpbegin, nproma);
                if constexpr (nproma_extent != dynamic_extent) {
                    constexpr idx_t nrof_static = nproma_extent;
                    copy_nonblocked_to_blocked_contiguous_rank3_block<nrof_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nproma, loop_order);
                }
                else {
                    copy_nonblocked_to_blocked_contiguous_rank3_block(
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
            if constexpr (nproma_extent != dynamic_extent) {
                constexpr idx_t nrof_static = nproma_extent;
                if (use_memcpy) {
                    std::memcpy(
                        raw_blocked_jblk,
                        raw_nonblocked_jblk,
                        static_cast<std::size_t>(nrof_static) * sizeof(Value));
                }
                else {
                    for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
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
void copy_blocked_to_nonblocked_contiguous_rank4_block(const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
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
                    ATLAS_RELAYOUT_SIMD
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        raw_nonblocked_jblk_jvar_jlev[jrof * point_stride] = raw_blocked_jblk_jvar_jlev[jrof];
                    }
                }
                else {
                    ATLAS_RELAYOUT_SIMD
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
                    ATLAS_RELAYOUT_SIMD
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
                    ATLAS_RELAYOUT_SIMD
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        raw_nonblocked_jblk_jrof_jvar[jlev * lev_nonblocked_stride] = raw_blocked_jblk_jrof_jvar[jlev * lev_block_stride];
                    }
                }
            }
        }
    }
}

template <idx_t nrof_static = 0, typename BlockedValue, typename NonblockedValue>
void copy_blocked_to_nonblocked_contiguous_rank3_block(const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
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
                ATLAS_RELAYOUT_SIMD
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    raw_nonblocked_jblk_jlev[jrof * nlev] = raw_blocked_jblk_jlev[jrof];
                }
            }
            else {
                ATLAS_RELAYOUT_SIMD
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
                ATLAS_RELAYOUT_SIMD
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    raw_nonblocked_jblk_jrof[jlev] = raw_blocked_jblk_jrof[jlev * nproma];
                }
            }
        }
        else {
            for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                const BlockedValue* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                NonblockedValue* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * nlev;
                ATLAS_RELAYOUT_SIMD
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    raw_nonblocked_jblk_jrof[jlev] = raw_blocked_jblk_jrof[jlev * nproma];
                }
            }
        }
    }
}

template <size_t nproma_extent, class Blocked, class Nonblocked>
void host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma(const Blocked blocked, Nonblocked nonblocked) {
    // Requirements: blocked and nonblocked must be contiguous, and if nproma_static is specified, it must match the last dimension of blocked.
    ATLAS_ASSERT(is_block_contiguous(blocked));
    ATLAS_ASSERT(has_layout_right(nonblocked));

    const idx_t nproma = last_extent(blocked);
    if constexpr (nproma_extent != dynamic_extent) {
        ATLAS_ASSERT(nproma_extent == nproma);
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
                const idx_t nrof = nproma;
                if constexpr (nproma_extent != dynamic_extent) {
                    constexpr idx_t nrof_static = nproma_extent;
                    copy_blocked_to_nonblocked_contiguous_rank4_block<nrof_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nvar, nproma, loop_order);
                }
                else {
                    copy_blocked_to_nonblocked_contiguous_rank4_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nvar, nproma, loop_order);
                }
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                if constexpr (nproma_extent != dynamic_extent) {
                    copy_blocked_to_nonblocked_contiguous_rank4_block<nproma_extent>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nvar, nproma_extent, loop_order);
                }
                else {
                    copy_blocked_to_nonblocked_contiguous_rank4_block(
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
                const idx_t nrof = nproma;
                if constexpr (nproma_extent != dynamic_extent) {
                    constexpr idx_t nrof_static = nproma_extent;
                    copy_blocked_to_nonblocked_contiguous_rank3_block<nrof_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nproma, loop_order);
                }
                else {

                    copy_blocked_to_nonblocked_contiguous_rank3_block(
                        raw_blocked_jblk, raw_nonblocked_jblk, nproma, nlev, nproma, loop_order);
                }
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                if constexpr (nproma_extent != dynamic_extent) {
                    constexpr idx_t nrof_static = nproma_extent;
                    copy_blocked_to_nonblocked_contiguous_rank3_block<nrof_static>(
                        raw_blocked_jblk, raw_nonblocked_jblk, nrof, nlev, nproma, loop_order);
                }
                else {
                    copy_blocked_to_nonblocked_contiguous_rank3_block(
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

#endif

// Fallback implementation wrapper with static nproma dispatch
template <typename Blocked, size_t nproma, BlockAlignment block_alignment, idx_t Rank = Blocked::rank()>
struct BlockedSubspan;

template <typename Blocked, size_t nproma, BlockAlignment block_alignment>
struct BlockedSubspan<Blocked, nproma, block_alignment, 4> {
    static constexpr bool is_aligned = (block_alignment == BlockAlignment::aligned);
    using value_t = mdspan_introspection_detail::view_value_t<Blocked>;
    static constexpr std::size_t alignment = (not is_aligned) ? alignof(value_t) : ::atlas::alignment;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t, alignment>;
    using extents_t = extents<idx_t, dynamic_extent, dynamic_extent, nproma>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};

template <typename Blocked, size_t nproma, BlockAlignment block_alignment>
struct BlockedSubspan<Blocked, nproma, block_alignment, 3> {
    static constexpr bool is_aligned = (block_alignment == BlockAlignment::aligned);
    using value_t = mdspan_introspection_detail::view_value_t<Blocked>;
    static constexpr std::size_t alignment = (not is_aligned) ? alignof(value_t) : ::atlas::alignment;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t, alignment>;
    using extents_t = extents<idx_t, dynamic_extent, nproma>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};

template <typename Blocked, size_t nproma, BlockAlignment block_alignment>
struct BlockedSubspan<Blocked, nproma, block_alignment, 2> {
    static constexpr bool is_aligned = (block_alignment == BlockAlignment::aligned);
    using value_t = mdspan_introspection_detail::view_value_t<Blocked>;
    static constexpr std::size_t alignment = (nproma == dynamic_extent || not is_aligned) ? alignof(value_t) : blocked_subspan_alignment_v<nproma, value_t>;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t, alignment>;
    using extents_t = extents<idx_t, nproma>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};

template <typename Nonblocked, size_t nproma, idx_t Rank = Nonblocked::rank()>
struct NonblockedSubspan;

template <typename Nonblocked, size_t nproma>
struct NonblockedSubspan<Nonblocked, nproma, 3> {
    using value_t = mdspan_introspection_detail::view_value_t<Nonblocked>;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t,64>;
    using extents_t = extents<idx_t, nproma, dynamic_extent, dynamic_extent>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};

template <typename Nonblocked, size_t nproma>
struct NonblockedSubspan<Nonblocked, nproma, 2> {
    using value_t = mdspan_introspection_detail::view_value_t<Nonblocked>;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t,64>;
    using extents_t = extents<idx_t, nproma, dynamic_extent>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};

template <typename Nonblocked, size_t nproma>
struct NonblockedSubspan<Nonblocked, nproma, 1> {
    using value_t = mdspan_introspection_detail::view_value_t<Nonblocked>;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t,64>;
    using extents_t = extents<idx_t, nproma>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};


template <typename Blocked, size_t nproma_extent = dynamic_extent, BlockAlignment block_alignment = BlockAlignment::aligned>
using blocked_subspan_t = typename BlockedSubspan<Blocked,nproma_extent,block_alignment>::span_t;

template <typename Nonblocked, size_t nproma = dynamic_extent>
using nonblocked_subspan_t = typename NonblockedSubspan<Nonblocked,nproma>::span_t;

template <size_t nproma = dynamic_extent, typename Blocked>
auto make_blocked_subspan_extents(const Blocked& blocked) {
    using extents_type = typename blocked_subspan_t<Blocked,nproma>::extents_type;
    if constexpr (Blocked::rank() == 4) {
        if constexpr (nproma == dynamic_extent) {
            return extents_type{blocked.extent(1), blocked.extent(2), blocked.extent(3)};
        }
        else {
            return extents_type{blocked.extent(1), blocked.extent(2)};
        }
    }
    else if constexpr (Blocked::rank() == 3) {
        if constexpr (nproma == dynamic_extent) {
            return extents_type{blocked.extent(1), blocked.extent(2)};
        }
        else {
            return extents_type{blocked.extent(1)};
        }
    }
    else if constexpr (Blocked::rank() == 2) {
        if constexpr (nproma == dynamic_extent) {
            return extents_type{blocked.extent(1)};
        }
        else {
            return extents_type{};
        }
    }
    else {
        throw_Exception("make_blocked_subspan_extents: Blocked rank not supported");
    }
}

template <size_t nproma = dynamic_extent, typename Nonblocked>
auto make_nonblocked_subspan_extents(const Nonblocked& nonblocked) {
    using extents_type = typename nonblocked_subspan_t<Nonblocked,nproma>::extents_type;
    if constexpr (Nonblocked::rank() == 3) {
        if constexpr (nproma == dynamic_extent) {
            return extents_type{nonblocked.extent(0), nonblocked.extent(1), nonblocked.extent(2)};
        }
        else {
            return extents_type{nonblocked.extent(1), nonblocked.extent(2)};
        }
    }
    else if constexpr (Nonblocked::rank() == 2) {
        if constexpr (nproma == dynamic_extent) {
            return extents_type{nonblocked.extent(0), nonblocked.extent(1)};
        }
        else {
            return extents_type{nonblocked.extent(1)};
        }
    }
    else if constexpr (Nonblocked::rank() == 1) {
        if constexpr (nproma == dynamic_extent) {
            return extents_type{nonblocked.extent(0)};
        }
        else {
            return extents_type{};
        }
    }
    else {
        throw_Exception("make_nonblocked_subspan_extents: Nonblocked rank not supported");
    }
}

template<typename Blocked>
bool is_block_aligned(const Blocked& block) {
    using Value = mdspan_introspection_detail::view_value_t<decltype(block)>;
    return is_aligned(block,alignment) &&
           (static_cast<std::size_t>(block.stride(1)) * sizeof(Value) % alignment == 0);
}

template<typename Blocked>
void assert_requirements_on_blocked(const Blocked& blocked) {
    ATLAS_ASSERT(is_block_contiguous(blocked));
}
template<typename Nonblocked>
void assert_requirements_on_nonblocked(const Nonblocked& nonblocked) {
    bool nonblocked_is_aligned = is_aligned(nonblocked,alignment);
    ATLAS_ASSERT(nonblocked_is_aligned);
    ATLAS_ASSERT(is_contiguous(nonblocked));
}

template <size_t nproma_extent, BlockAlignment block_alignment, class Blocked, class Nonblocked>
struct CopyBlockedToNonblockedBlock {
    using blocked_subspan_type = blocked_subspan_t<Blocked, nproma_extent, block_alignment>;
    using nonblocked_subspan_type = nonblocked_subspan_t<Nonblocked, nproma_extent>;
    using block_extents_t = typename blocked_subspan_type::extents_type;
    using nonblocked_extents_t = typename nonblocked_subspan_type::extents_type;

    static constexpr idx_t static_nrof() {
        if constexpr(nproma_extent != dynamic_extent) {
            return nproma_extent;
        }
        else {
            return 0;
        }
    }

    Blocked blocked;
    mutable Nonblocked nonblocked;
    idx_t np;
    idx_t nblks;
    idx_t nproma;
    RelayoutLoopOrder loop_order;
    block_extents_t block_extents;
    nonblocked_extents_t nonblocked_extents;

    CopyBlockedToNonblockedBlock(Blocked blocked, Nonblocked nonblocked):
        blocked(blocked),
        nonblocked(nonblocked),
        np(nonblocked.extent(0)),
        nblks(blocked.extent(0)),
        nproma(last_extent(blocked)),
        loop_order(relayout_loop_order()),
        block_extents(make_blocked_subspan_extents<nproma_extent>(blocked)),
        nonblocked_extents(make_nonblocked_subspan_extents<nproma_extent>(nonblocked)) {}

    void operator()(idx_t jblk) const {
        const idx_t jpbegin = jblk * nproma;

        if constexpr(Blocked::rank()==4) {
            blocked_subspan_type block_jblk{&blocked(jblk, 0, 0, 0), block_extents};
            nonblocked_subspan_type nonblocked_jblk{&nonblocked(jpbegin, 0, 0), nonblocked_extents};
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_blocked_to_nonblocked_rank4_block<static_nrof()>(block_jblk, nonblocked_jblk, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_blocked_to_nonblocked_rank4_block(block_jblk, nonblocked_jblk, nrof);
            }
        }
        else if constexpr(Blocked::rank()==3) {
            blocked_subspan_type block_jblk{&blocked(jblk, 0, 0), block_extents};
            nonblocked_subspan_type nonblocked_jblk{&nonblocked(jpbegin, 0), nonblocked_extents};
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_blocked_to_nonblocked_rank3_block<static_nrof()>(block_jblk, nonblocked_jblk, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_blocked_to_nonblocked_rank3_block(block_jblk, nonblocked_jblk, nrof);
            }
        }
        else if constexpr (Blocked::rank()==2) {
            blocked_subspan_type block_jblk{&blocked(jblk, 0), block_extents};
            nonblocked_subspan_type nonblocked_jblk{&nonblocked(jpbegin), nonblocked_extents};
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    nonblocked_jblk(jrof) = block_jblk(jrof);
                }
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    nonblocked_jblk(jrof) = block_jblk(jrof);
                }
            }
        }
        else {
            ATLAS_THROW_EXCEPTION("host_copy_blocked_to_nonblocked_mdspan not implemented");
        }
    }

private:

    template <idx_t nrof_static = 0>
    void copy_blocked_to_nonblocked_rank4_block(blocked_subspan_type& block,
                                                const nonblocked_subspan_type& nonblocked_chunk,
                                                [[maybe_unused]] const idx_t nrof) const {
        if constexpr (nrof_static == 0) {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jvar = 0; jvar < block.extent(0); ++jvar) {
                    for (idx_t jlev = 0; jlev < block.extent(1); ++jlev) {
                        for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                            nonblocked_chunk(jrof,jlev,jvar) = block(jvar,jlev,jrof);
                        }
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    for (idx_t jlev = 0; jlev < block.extent(1); ++jlev) {
                        for (idx_t jvar = 0; jvar < block.extent(0); ++jvar) {
                            nonblocked_chunk(jrof,jlev,jvar) = block(jvar,jlev,jrof);
                        }
                    }
                }
            }
        }
        else {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jvar = 0; jvar < block.extent(0); ++jvar) {
                    for (idx_t jlev = 0; jlev < block.extent(1); ++jlev) {
                        for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                            nonblocked_chunk(jrof,jlev,jvar) = block(jvar,jlev,jrof);
                        }
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    for (idx_t jlev = 0; jlev < block.extent(1); ++jlev) {
                        for (idx_t jvar = 0; jvar < block.extent(0); ++jvar) {
                            nonblocked_chunk(jrof,jlev,jvar) = block(jvar,jlev,jrof);
                        }
                    }
                }
            }
        }
    }

    template <idx_t nrof_static = 0>
    void copy_blocked_to_nonblocked_rank3_block(blocked_subspan_type& block,
                                                const nonblocked_subspan_type& nonblocked_chunk,
                                                [[maybe_unused]] const idx_t nrof) const {
        if constexpr (nrof_static == 0) {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jlev = 0; jlev < block.extent(0); ++jlev) {
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        nonblocked_chunk(jrof,jlev) = block(jlev,jrof);
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    for (idx_t jlev = 0; jlev < block.extent(0); ++jlev) {
                        nonblocked_chunk(jrof,jlev) = block(jlev,jrof);
                    }
                }
            }
        }
        else {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jlev = 0; jlev < block.extent(0); ++jlev) {
                    for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                        nonblocked_chunk(jrof,jlev) = block(jlev,jrof);
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    for (idx_t jlev = 0; jlev < block.extent(0); ++jlev) {
                        nonblocked_chunk(jrof,jlev) = block(jlev,jrof);
                    }
                }
            }
        }
    }
};

template <size_t nproma_extent, BlockAlignment block_alignment, class Blocked, class Nonblocked>
auto make_copy_blocked_to_nonblocked_block(const Blocked blocked, Nonblocked nonblocked) {
    return CopyBlockedToNonblockedBlock<nproma_extent, block_alignment, Blocked, Nonblocked>{blocked, nonblocked};
}

template <size_t nproma_extent, BlockAlignment block_alignment, class Nonblocked, class Blocked>
struct CopyNonblockedToBlockedBlock {
    using blocked_subspan_type = blocked_subspan_t<Blocked, nproma_extent, block_alignment>;
    using nonblocked_subspan_type = nonblocked_subspan_t<Nonblocked, nproma_extent>;
    using block_extents_t = typename blocked_subspan_type::extents_type;
    using nonblocked_extents_t = typename nonblocked_subspan_type::extents_type;

    static constexpr idx_t static_nrof() {
        if constexpr(nproma_extent != dynamic_extent) {
            return nproma_extent;
        }
        else {
            return 0;
        }
    }

    Nonblocked nonblocked;
    mutable Blocked blocked;
    idx_t np;
    idx_t nblks;
    idx_t nproma;
    RelayoutLoopOrder loop_order;
    block_extents_t block_extents;
    nonblocked_extents_t nonblocked_extents;

    CopyNonblockedToBlockedBlock(Nonblocked nonblocked, Blocked blocked):
        nonblocked(nonblocked),
        blocked(blocked),
        np(nonblocked.extent(0)),
        nblks(blocked.extent(0)),
        nproma(last_extent(blocked)),
        loop_order(relayout_loop_order()),
        block_extents(make_blocked_subspan_extents<nproma_extent>(blocked)),
        nonblocked_extents(make_nonblocked_subspan_extents<nproma_extent>(nonblocked)) {}

    void operator()(idx_t jblk) const {
        const idx_t jpbegin = jblk * nproma;

        if constexpr(Blocked::rank()==4) {
            blocked_subspan_type block_jblk{&blocked(jblk, 0, 0, 0), block_extents};
            nonblocked_subspan_type nonblocked_jblk{&nonblocked(jpbegin, 0, 0), nonblocked_extents};
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_nonblocked_to_blocked_rank4_block<static_nrof()>(nonblocked_jblk, block_jblk, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_nonblocked_to_blocked_rank4_block(nonblocked_jblk, block_jblk, nrof);
            }
        }
        else if constexpr(Blocked::rank()==3) {
            blocked_subspan_type block_jblk{&blocked(jblk, 0, 0), block_extents};
            nonblocked_subspan_type nonblocked_jblk{&nonblocked(jpbegin, 0), nonblocked_extents};
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_nonblocked_to_blocked_rank3_block<static_nrof()>(nonblocked_jblk, block_jblk, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_nonblocked_to_blocked_rank3_block(nonblocked_jblk, block_jblk, nrof);
            }
        }
        else if constexpr (Blocked::rank()==2) {
            blocked_subspan_type block_jblk{&blocked(jblk, 0), block_extents};
            nonblocked_subspan_type nonblocked_jblk{&nonblocked(jpbegin), nonblocked_extents};
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    block_jblk(jrof) = nonblocked_jblk(jrof);
                }
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    block_jblk(jrof) = nonblocked_jblk(jrof);
                }
            }
        }
        else {
            ATLAS_THROW_EXCEPTION("host_copy_nonblocked_to_blocked_mdspan not implemented");
        }
    }

private:

    template <idx_t nrof_static = 0>
    void copy_nonblocked_to_blocked_rank4_block(const nonblocked_subspan_type& nonblocked_chunk, blocked_subspan_type& block,
                                                [[maybe_unused]] const idx_t nrof) const {
        if constexpr (nrof_static == 0) {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jvar = 0; jvar < block.extent(0); ++jvar) {
                    for (idx_t jlev = 0; jlev < block.extent(1); ++jlev) {
                        for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                            block(jvar,jlev,jrof) = nonblocked_chunk(jrof,jlev,jvar);
                        }
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    for (idx_t jlev = 0; jlev < block.extent(1); ++jlev) {
                        for (idx_t jvar = 0; jvar < block.extent(0); ++jvar) {
                            block(jvar,jlev,jrof) = nonblocked_chunk(jrof,jlev,jvar);
                        }
                    }
                }
            }
        }
        else {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jvar = 0; jvar < block.extent(0); ++jvar) {
                    for (idx_t jlev = 0; jlev < block.extent(1); ++jlev) {
                        for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                            block(jvar,jlev,jrof) = nonblocked_chunk(jrof,jlev,jvar);
                        }
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    for (idx_t jlev = 0; jlev < block.extent(1); ++jlev) {
                        for (idx_t jvar = 0; jvar < block.extent(0); ++jvar) {
                            block(jvar,jlev,jrof) = nonblocked_chunk(jrof,jlev,jvar);
                        }
                    }
                }
            }
        }
    }

    template <idx_t nrof_static = 0>
    void copy_nonblocked_to_blocked_rank3_block(const nonblocked_subspan_type& nonblocked_chunk, blocked_subspan_type& block,
                                                [[maybe_unused]] const idx_t nrof) const {
        if constexpr (nrof_static == 0) {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jlev = 0; jlev < block.extent(0); ++jlev) {
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        block(jlev,jrof) = nonblocked_chunk(jrof,jlev);
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    for (idx_t jlev = 0; jlev < block.extent(0); ++jlev) {
                        block(jlev,jrof) = nonblocked_chunk(jrof,jlev);
                    }
                }
            }
        }
        else {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jlev = 0; jlev < block.extent(0); ++jlev) {
                    for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                        block(jlev,jrof) = nonblocked_chunk(jrof,jlev);
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    for (idx_t jlev = 0; jlev < block.extent(0); ++jlev) {
                        block(jlev,jrof) = nonblocked_chunk(jrof,jlev);
                    }
                }
            }
        }
    }

};

template <size_t nproma_extent, BlockAlignment block_alignment, class Nonblocked, class Blocked>
auto make_copy_nonblocked_to_blocked_block(const Nonblocked nonblocked, Blocked blocked) {
    return CopyNonblockedToBlockedBlock<nproma_extent, block_alignment, Nonblocked, Blocked>{nonblocked, blocked};
}

template <size_t nproma_extent, class Blocked, class Nonblocked>
void host_copy_blocked_to_nonblocked_nproma(const Blocked blocked, Nonblocked nonblocked) {
    assert_requirements_on_blocked(blocked);
    assert_requirements_on_nonblocked(nonblocked);

    if (is_block_aligned(blocked)) {
        auto copy_blocked_to_nonblocked_block = make_copy_blocked_to_nonblocked_block<nproma_extent, BlockAlignment::aligned>(blocked, nonblocked);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < copy_blocked_to_nonblocked_block.nblks; ++jblk) {
            copy_blocked_to_nonblocked_block(jblk);
        }
    }
    else {
        Log::debug() << "host_copy_blocked_to_nonblocked_nproma: A block is not aligned, falling back to unaligned copy.";
        Log::debug() << "\nBlocked: " << std::vector<idx_t>(blocked.shape(), blocked.shape() + blocked.rank()) << std::endl;

        auto copy_blocked_to_nonblocked_block = make_copy_blocked_to_nonblocked_block<nproma_extent, BlockAlignment::unaligned>(blocked, nonblocked);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < copy_blocked_to_nonblocked_block.nblks; ++jblk) {
            copy_blocked_to_nonblocked_block(jblk);
        }
    }
}

template <size_t nproma_extent, class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_nproma(const Nonblocked nonblocked, Blocked blocked) {
    assert_requirements_on_nonblocked(nonblocked);
    assert_requirements_on_blocked(blocked);

    if (is_block_aligned(blocked)) {
        auto copy_nonblocked_to_blocked_block = make_copy_nonblocked_to_blocked_block<nproma_extent, BlockAlignment::aligned>(nonblocked, blocked);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < copy_nonblocked_to_blocked_block.nblks; ++jblk) {
            copy_nonblocked_to_blocked_block(jblk);
        }
    }
    else {
        Log::debug() << "host_copy_nonblocked_to_blocked_nproma: A block is not aligned, falling back to unaligned copy.";
        Log::debug() << "\nBlocked: " << std::vector<idx_t>(blocked.shape(), blocked.shape() + blocked.rank()) << std::endl;

        auto copy_nonblocked_to_blocked_block = make_copy_nonblocked_to_blocked_block<nproma_extent, BlockAlignment::unaligned>(nonblocked, blocked);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < copy_nonblocked_to_blocked_block.nblks; ++jblk) {
            copy_nonblocked_to_blocked_block(jblk);
        }
    }
}

}  // namespace

// Fallback implementation wrappers with static nproma dispatch
template <idx_t nproma, class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_fallback_nproma(const Nonblocked nonblocked, Blocked blocked) {
    const idx_t np    = nonblocked.extent(0);
    const idx_t nblks = blocked.extent(0);
    [[maybe_unused]] const RelayoutLoopOrder loop_order = relayout_loop_order();

    if constexpr(blocked.rank()==4) {
        const idx_t nlev = nonblocked.extent(1);
        const idx_t nvar = nonblocked.extent(2);

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
        const idx_t nlev = nonblocked.extent(1);

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
        ATLAS_THROW_EXCEPTION("host_copy_nonblocked_to_blocked_mdspan not implemented");
    }
}

template <class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_fallback(const Nonblocked nonblocked, Blocked blocked) {
    const idx_t np    = nonblocked.extent(0);
    const idx_t nblks = blocked.extent(0);
    const idx_t nproma = blocked.extent(blocked.rank()-1);
    [[maybe_unused]] const RelayoutLoopOrder loop_order = relayout_loop_order();

    if constexpr(blocked.rank()==4) {
        const idx_t nlev = nonblocked.extent(1);
        const idx_t nvar = nonblocked.extent(2);

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
        const idx_t nlev = nonblocked.extent(1);

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
        ATLAS_THROW_EXCEPTION("host_copy_nonblocked_to_blocked_mdspan not implemented");
    }
}

template <class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_impl(const Nonblocked nonblocked, Blocked blocked) {
    static_assert(nonblocked.rank() == blocked.rank()-1);
    const idx_t nproma = blocked.extent(blocked.rank() - 1);

    // If Blocked is an mdspan whose last (nproma) dimension is a static extent, then nproma is
    // known at compile time.  In that case we can call the statically-sized kernels directly and
    // skip the runtime switch dispatch below, avoiding the unused template instantiations.
    constexpr std::size_t nproma_extent = []() {
        if constexpr (is_mdspan<Blocked>::value) {
            return Blocked::static_extent(Blocked::rank() - 1);
        }
        else {
            return dynamic_extent;
        }
    }();

#if DISABLE_RAW_POINTERS == 0
    if (relayout_implementation() == RelayoutImplementation::raw_pointers && has_layout_right(nonblocked) && is_block_contiguous(blocked)) {
        if (relayout_nproma_dispatch() != RelayoutNpromaDispatch::static_dispatch) {
            return host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma<dynamic_extent>(nonblocked, blocked);
        }

        // Optimized paths possible with static nproma dispatch.
        if constexpr (nproma_extent != dynamic_extent) {
            return host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma<nproma_extent>(nonblocked, blocked);
        }
        else {
            switch (nproma) {
                case 8:   return host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma<8 >(nonblocked, blocked);
                case 16:  return host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma<16>(nonblocked, blocked);
                case 32:  return host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma<32>(nonblocked, blocked);
                case 64:  return host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma<64>(nonblocked, blocked);
                case 128: return host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma<128>(nonblocked, blocked);
                case 256: return host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma<256>(nonblocked, blocked);
                default:  return host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma<dynamic_extent>(nonblocked, blocked);
            }
        }
    }
#endif

    // Fall back to generic implementation if the views are not contiguous or block-contiguous.
    // Apply static dispatch even for fallback paths for better performance.
    if (relayout_nproma_dispatch() != RelayoutNpromaDispatch::static_dispatch) {
        // Runtime dispatch for fallback
        host_copy_nonblocked_to_blocked_nproma<dynamic_extent>(nonblocked, blocked);
    }

    // Optimized fallback paths with static nproma dispatch.
    if constexpr (nproma_extent != dynamic_extent) {
        return host_copy_nonblocked_to_blocked_nproma<nproma_extent>(nonblocked, blocked);
    }
    else {
        switch (nproma) {
            case 8:   return host_copy_nonblocked_to_blocked_nproma<8  >(nonblocked, blocked);
            case 16:  return host_copy_nonblocked_to_blocked_nproma<16 >(nonblocked, blocked);
            case 32:  return host_copy_nonblocked_to_blocked_nproma<32 >(nonblocked, blocked);
            case 64:  return host_copy_nonblocked_to_blocked_nproma<64 >(nonblocked, blocked);
            case 128: return host_copy_nonblocked_to_blocked_nproma<128>(nonblocked, blocked);
            case 256: return host_copy_nonblocked_to_blocked_nproma<256>(nonblocked, blocked);
            default:  return host_copy_nonblocked_to_blocked_nproma<dynamic_extent>(nonblocked, blocked);
        }
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
void host_copy_blocked_to_nonblocked_impl(const Blocked blocked, Nonblocked nonblocked) {
    static_assert(nonblocked.rank() == blocked.rank()-1);

    // If Blocked is an mdspan whose last (nproma) dimension is a static extent, then nproma is
    // known at compile time.  In that case we can call the statically-sized kernels directly and
    // skip the runtime switch dispatch below, avoiding the unused template instantiations.
    constexpr std::size_t nproma_extent = []() {
        if constexpr (is_mdspan<Blocked>::value) {
            return Blocked::static_extent(Blocked::rank() - 1);
        }
        else {
            return dynamic_extent;
        }
    }();

#if DISABLE_RAW_POINTERS == 0
    if (relayout_implementation() == RelayoutImplementation::raw_pointers && has_layout_right(nonblocked) && is_block_contiguous(blocked)) {
        if (relayout_nproma_dispatch() != RelayoutNpromaDispatch::static_dispatch) {
            return host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma<dynamic_extent>(blocked, nonblocked);
        }

        // Optimized paths possible with static nproma dispatch.
        if constexpr (nproma_extent != dynamic_extent) {
            return host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma<nproma_extent>(blocked, nonblocked);
        }
        else {
            auto nproma = last_extent(blocked);
            switch (nproma) {
                case 8:   return host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma<8  >(blocked, nonblocked);
                case 16:  return host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma<16 >(blocked, nonblocked);
                case 32:  return host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma<32 >(blocked, nonblocked);
                case 64:  return host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma<64 >(blocked, nonblocked);
                case 128: return host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma<128>(blocked, nonblocked);
                case 256: return host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma<256>(blocked, nonblocked);
                default:  return host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma<dynamic_extent>(blocked, nonblocked);
            }
        }
    }
#endif

    // Fall back to generic implementation if the views are not contiguous or block-contiguous.
    // Apply static dispatch even for fallback paths for better performance.
    if (relayout_nproma_dispatch() != RelayoutNpromaDispatch::static_dispatch) {
        // Runtime dispatch for fallback
        return host_copy_blocked_to_nonblocked_nproma<dynamic_extent>(blocked, nonblocked);
    }

    // Optimized fallback paths with static nproma dispatch.
    if constexpr (nproma_extent != dynamic_extent) {
        return host_copy_blocked_to_nonblocked_nproma<nproma_extent>(blocked, nonblocked);
    }
    else {
        auto nproma = last_extent(blocked);
        switch (nproma) {
            case 8:   return host_copy_blocked_to_nonblocked_nproma<8  >(blocked, nonblocked);
            case 16:  return host_copy_blocked_to_nonblocked_nproma<16 >(blocked, nonblocked);
            case 32:  return host_copy_blocked_to_nonblocked_nproma<32 >(blocked, nonblocked);
            case 64:  return host_copy_blocked_to_nonblocked_nproma<64 >(blocked, nonblocked);
            case 128: return host_copy_blocked_to_nonblocked_nproma<128>(blocked, nonblocked);
            case 256: return host_copy_blocked_to_nonblocked_nproma<256>(blocked, nonblocked);
            default:  return host_copy_blocked_to_nonblocked_nproma<dynamic_extent>(blocked, nonblocked);
        }
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
void host_copy_blocked_to_blocked_impl(const BlockedIn blocked_in, BlockedOut blocked_out) {
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

template <class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked) {
    return host_copy_nonblocked_to_blocked_impl(nonblocked, blocked);
}

template <class Blocked, class Nonblocked>
void host_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    return host_copy_blocked_to_nonblocked_impl(blocked, nonblocked);
}

template <class BlockedIn, class BlockedOut>
void host_copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out) {
    return host_copy_blocked_to_blocked_impl(blocked_in, blocked_out);
}

#undef ATLAS_RELAYOUT_RESTRICT

}  // namespace atlas


#define ATLAS_RELAYOUT_EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(TYPE, BLOCKED_RANK) \
    template void atlas::host_copy_blocked_to_nonblocked_mdspan<atlas::array::ArrayView<const TYPE,BLOCKED_RANK>,atlas::array::ArrayView<TYPE,BLOCKED_RANK-1>>(atlas::array::ArrayView<const TYPE,BLOCKED_RANK>, atlas::array::ArrayView<TYPE,BLOCKED_RANK-1>); \
    template void atlas::host_copy_nonblocked_to_blocked_mdspan<atlas::array::ArrayView<const TYPE,BLOCKED_RANK-1>,atlas::array::ArrayView<TYPE,BLOCKED_RANK>>(atlas::array::ArrayView<const TYPE,BLOCKED_RANK-1>, atlas::array::ArrayView<TYPE,BLOCKED_RANK>); \
    template void atlas::host_copy_blocked_to_blocked_mdspan<atlas::array::ArrayView<const TYPE,BLOCKED_RANK>,atlas::array::ArrayView<TYPE,BLOCKED_RANK>>(atlas::array::ArrayView<const TYPE,BLOCKED_RANK>, atlas::array::ArrayView<TYPE,BLOCKED_RANK>); \
