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
#ifndef ATLAS_RELAYOUT_SIMD
#define ATLAS_RELAYOUT_SIMD
#endif

// Independently toggleable SIMD directive for the mdspan implementation kernels. The strided
// (transpose) inner loops of the mdspan rank-3/rank-4 kernels do not auto-vectorize without an
// explicit directive; define this to `atlas_omp_pragma(omp simd)` to force vectorization, or
// comment it out to leave those loops scalar.
#define ATLAS_RELAYOUT_MDSPAN_SIMD atlas_omp_pragma(omp simd)
#ifndef ATLAS_RELAYOUT_MDSPAN_SIMD
#define ATLAS_RELAYOUT_MDSPAN_SIMD
#endif

#define ATLAS_RELAYOUT_PROFILING 0
#if ATLAS_RELAYOUT_PROFILING
#define ATLAS_RELAYOUT_NOINLINE_IF_PROFILING_IF_PROFILING __attribute__((noinline))
#else
#define ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
#endif

#include "atlas/util/relayout.h"

#include <new>
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

#include "atlas/runtime/Trace.h"

using namespace atlas::array::introspection;

#ifdef atlas_omp_parallel_for
#undef atlas_omp_parallel_for
#endif
#define atlas_omp_parallel_for for
namespace atlas {

namespace {

// SIMD/vectorization alignment (in bytes) used for aligned accessor policies and for aligning
// blocked storage so that inner `nproma` loops can use aligned vector loads/stores.  This is the
// widest vector-register width the target can access aligned -- NOT the cache-line size.  We select
// it from the compiler's SIMD feature macros, and fall back to 16 bytes (SSE2 / ARM NEON, including
// Apple Silicon), which is the common minimum for contemporary SIMD.  Over-aligning is only a minor
// waste; under-aligning is caught at runtime by is_aligned() which falls back to unaligned kernels.
#if defined(__AVX512F__)
constexpr std::size_t alignment = 64;   // AVX-512: 512-bit vectors
#elif defined(__AVX__)
constexpr std::size_t alignment = 32;   // AVX / AVX2: 256-bit vectors
#else
constexpr std::size_t alignment = 16;   // SSE2 and ARM NEON (incl. Apple Silicon): 128-bit vectors
#endif

constexpr std::size_t greatest_common_divisor(std::size_t lhs, std::size_t rhs) {
    while (rhs != 0) {
        const std::size_t remainder = lhs % rhs;
        lhs = rhs;
        rhs = remainder;
    }
    return lhs;
}

template <size_t nproma, typename ElementType>
inline constexpr std::size_t blocked_subspan_alignment_v =
    greatest_common_divisor(alignment, static_cast<std::size_t>(nproma) * sizeof(ElementType));

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

// TODO(willem): Add support for a benchmark-generated host relayout tuning profile
// that can be read once at runtime and used here to select implementation and/or
// loop order for shape buckets on the current system, while still letting explicit
// environment settings take precedence.

[[maybe_unused]] RelayoutLoopOrder relayout_loop_order() {
    static RelayoutLoopOrder cached_loop_order = []() {
         const char* loop_order = std::getenv(relayout_loop_order_env);
         if (loop_order && std::strcmp(loop_order, "nproma_outermost") == 0) {
             return RelayoutLoopOrder::nproma_outermost;
         }
         // nonblocked_coalesced is a device-only variant; on the host it maps to nproma_outermost.
         if (loop_order && std::strcmp(loop_order, "nonblocked_coalesced") == 0) {
             return RelayoutLoopOrder::nproma_outermost;
         }
         // coalesced_write / coalesced_read are device-only meta-orders; on the host, where there
         // is no memory coalescing, they map to nproma_outermost.
         if (loop_order && (std::strcmp(loop_order, "coalesced_write") == 0 ||
                            std::strcmp(loop_order, "coalesced_read") == 0)) {
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

template <typename View, typename = void>
struct get_nproma_extent {
    static constexpr std::size_t value = dynamic_extent;
};

template <typename View>
struct get_nproma_extent<View, std::enable_if_t<array::introspection::is_mdspan<View>()>> {
    static constexpr std::size_t value = View::static_extent(View::rank() - 1);
};

template <typename View>
static inline constexpr std::size_t nproma_extent_v = get_nproma_extent<View>::value;



#if DISABLE_RAW_POINTERS == 0
template <size_t nproma_extent, class Nonblocked, class Blocked>
/// @brief Per-block raw-pointer kernel for nonblocked-to-blocked host relayout.
///
/// This worker is selected when the nonblocked input is layout-right and the blocked
/// output is block-contiguous, allowing direct pointer arithmetic rather than mdspan
/// subviews. `nproma_extent` may be a compile-time block width or `dynamic_extent`
/// when the tail width is known only at runtime.
///
/// `operator()(jblk)` copies exactly one blocked chunk, using the full `nproma`
/// width for interior blocks and a shortened width for the final partial block.
struct CopyNonblockedToBlockedContiguousRawPointers {
    using blocked_element_type = array::introspection::element_t<Blocked>;
    using nonblocked_element_type = array::introspection::element_t<Nonblocked>;
    using value_type = std::remove_cv_t<blocked_element_type>;

    static constexpr idx_t static_nrof() {
        if constexpr (nproma_extent != dynamic_extent) {
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
    bool use_memcpy;

    CopyNonblockedToBlockedContiguousRawPointers(Nonblocked nonblocked, Blocked blocked):
        nonblocked(nonblocked),
        blocked(blocked),
        np(nonblocked.extent(0)),
        nblks(blocked.extent(0)),
        nproma(last_extent(blocked)),
        loop_order(relayout_loop_order()),
        use_memcpy(relayout_blocked_nonblocked_use_memcpy()) {
        if constexpr (Blocked::rank() == 4) {
            ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
            ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
        }
        else if constexpr (Blocked::rank() == 3) {
            ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
        }
    }

    /// @brief Copy one logical blocked chunk into the output blocked view.
    /// @param jblk Block index in the blocked output view.
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING 
    void operator()(idx_t jblk) const {
        const idx_t jpbegin = jblk * nproma;

        if constexpr(Blocked::rank()==4) {
            auto* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk = &blocked(jblk, 0, 0, 0);
            const auto* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk = &nonblocked(jpbegin, 0, 0);
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_nonblocked_to_blocked_rank4_block<static_nrof()>(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_nonblocked_to_blocked_rank4_block(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
        }
        else if constexpr (Blocked::rank()==3) {
            auto* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk = &blocked(jblk, 0, 0);
            const auto* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk = &nonblocked(jpbegin, 0);
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_nonblocked_to_blocked_rank3_block<static_nrof()>(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_nonblocked_to_blocked_rank3_block(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
        }
        else if constexpr (Blocked::rank()==2) {
            auto* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk = &blocked(jblk, 0);
            const auto* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk = &nonblocked(jpbegin);
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_nonblocked_to_blocked_rank2_block<static_nrof()>(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_nonblocked_to_blocked_rank2_block(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
        }
        else {
            ATLAS_THROW_EXCEPTION("host_copy_nonblocked_to_blocked_mdspan not implemented");
        }
    }

private:
    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_nonblocked_to_blocked_rank4_block(value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                                const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                                [[maybe_unused]] const idx_t nrof) const {
        const idx_t nlev = nonblocked.extent(1);
        const idx_t nvar = nonblocked.extent(2);
        const idx_t point_stride = nlev * nvar;
        const idx_t var_block_stride = nlev * nproma;
        const idx_t lev_block_stride = nproma;
        const idx_t lev_nonblocked_stride = nvar;

        if (loop_order == RelayoutLoopOrder::nproma_innermost) {
            for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jvar = raw_blocked_jblk + jvar * var_block_stride;
                const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jvar = raw_nonblocked_jblk + jvar;
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jvar_jlev = raw_blocked_jblk_jvar + jlev * lev_block_stride;
                    const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jvar_jlev = raw_nonblocked_jblk_jvar + jlev * lev_nonblocked_stride;
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
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                    const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * point_stride;
                    for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                        value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof_jvar = raw_blocked_jblk_jrof + jvar * var_block_stride;
                        const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof_jvar = raw_nonblocked_jblk_jrof + jvar;
                        ATLAS_RELAYOUT_SIMD
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            raw_blocked_jblk_jrof_jvar[jlev * lev_block_stride] = raw_nonblocked_jblk_jrof_jvar[jlev * lev_nonblocked_stride];
                        }
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                    const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * point_stride;
                    for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                        value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof_jvar = raw_blocked_jblk_jrof + jvar * var_block_stride;
                        const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof_jvar = raw_nonblocked_jblk_jrof + jvar;
                        ATLAS_RELAYOUT_SIMD
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            raw_blocked_jblk_jrof_jvar[jlev * lev_block_stride] = raw_nonblocked_jblk_jrof_jvar[jlev * lev_nonblocked_stride];
                        }
                    }
                }
            }
        }
    }

    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_nonblocked_to_blocked_rank3_block(value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                                const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                                [[maybe_unused]] const idx_t nrof) const {
        constexpr auto nproma_static = nrof_static;
        const idx_t nlev = nonblocked.extent(1);

        if (loop_order == RelayoutLoopOrder::nproma_innermost) {
            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jlev = raw_nonblocked_jblk + jlev;
                if constexpr (nrof_static == 0) {
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jlev = raw_blocked_jblk + jlev * nproma;
                    ATLAS_RELAYOUT_SIMD
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        raw_blocked_jblk_jlev[jrof] = raw_nonblocked_jblk_jlev[jrof * nlev];
                    }
                }
                else {
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jlev = raw_blocked_jblk + jlev * nproma_static;
                    ATLAS_RELAYOUT_SIMD
                    for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                        raw_blocked_jblk_jlev[jrof] = raw_nonblocked_jblk_jlev[jrof * nlev];
                    }
                }
            }
        }
        else {
            if constexpr (nrof_static == 0) {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                    const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * nlev;
                    ATLAS_RELAYOUT_SIMD
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        raw_blocked_jblk_jrof[jlev * nproma] = raw_nonblocked_jblk_jrof[jlev];
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                    const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * nlev;
                    ATLAS_RELAYOUT_SIMD
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        raw_blocked_jblk_jrof[jlev * nproma_static] = raw_nonblocked_jblk_jrof[jlev];
                    }
                }
            }
        }
    }

    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_nonblocked_to_blocked_rank2_block(value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                                const value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                                [[maybe_unused]] const idx_t nrof) const {
        if (use_memcpy) {
            const std::size_t count = static_cast<std::size_t>((nrof_static == 0) ? nrof : nrof_static);
            std::memcpy(raw_blocked_jblk, raw_nonblocked_jblk, count * sizeof(value_type));
        }
        else if constexpr (nrof_static == 0) {
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                raw_blocked_jblk[jrof] = raw_nonblocked_jblk[jrof];
            }
        }
        else {
            for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                raw_blocked_jblk[jrof] = raw_nonblocked_jblk[jrof];
            }
        }
    }
};

template <size_t nproma_extent, class Nonblocked, class Blocked>
auto make_copy_nonblocked_to_blocked_contiguous_raw_pointers(const Nonblocked nonblocked, Blocked blocked) {
    return CopyNonblockedToBlockedContiguousRawPointers<nproma_extent, Nonblocked, Blocked>{nonblocked, blocked};
}

template <size_t nproma_extent, class Nonblocked, class Blocked>
/// @brief Launch the raw-pointer nonblocked-to-blocked kernel for all blocks.
///
/// Preconditions:
/// - `nonblocked.rank() == blocked.rank() - 1`
/// - `blocked` is block-contiguous
/// - `nonblocked` behaves like layout-right contiguous storage
ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
void host_copy_nonblocked_to_blocked_contiguous_raw_pointers_nproma(const Nonblocked nonblocked, Blocked blocked) {
    ATLAS_ASSERT(is_block_contiguous(blocked));
    ATLAS_ASSERT(array::introspection::can_use_layout_right(nonblocked));

    const idx_t nblks  = blocked.extent(0);
    const idx_t nproma = last_extent(blocked);
    if constexpr (nproma_extent != dynamic_extent) {
        ATLAS_ASSERT(nproma_extent == nproma);
    }

    static_assert(nonblocked.rank() == blocked.rank()-1);

    auto copy_nonblocked_to_blocked_block = make_copy_nonblocked_to_blocked_contiguous_raw_pointers<nproma_extent>(nonblocked, blocked);
    atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
        copy_nonblocked_to_blocked_block(jblk);
    }
}

template <size_t nproma_extent, class Blocked, class Nonblocked>
/// @brief Per-block raw-pointer kernel for blocked-to-nonblocked host relayout.
///
/// This is the inverse of `CopyNonblockedToBlockedContiguousRawPointers`. It is used
/// when the blocked input is block-contiguous and the nonblocked output is layout-right,
/// so the copy can be expressed in terms of restricted raw pointers and simple strides.
struct CopyBlockedToNonblockedContiguousRawPointers {
    using blocked_element_type = array::introspection::element_t<Blocked>;
    using nonblocked_element_type = array::introspection::element_t<Nonblocked>;
    using value_type = std::remove_cv_t<nonblocked_element_type>;

    static constexpr idx_t static_nrof() {
        if constexpr (nproma_extent != dynamic_extent) {
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
    bool use_memcpy;

    CopyBlockedToNonblockedContiguousRawPointers(Blocked blocked, Nonblocked nonblocked):
        blocked(blocked),
        nonblocked(nonblocked),
        np(nonblocked.extent(0)),
        nblks(blocked.extent(0)),
        nproma(last_extent(blocked)),
        loop_order(relayout_loop_order()),
        use_memcpy(relayout_blocked_nonblocked_use_memcpy()) {
        if constexpr (Blocked::rank() == 4) {
            ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
            ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
        }
        else if constexpr (Blocked::rank() == 3) {
            ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
        }
    }

    /// @brief Copy one blocked chunk into its corresponding logical point range.
    /// @param jblk Block index in the blocked input view.
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void operator()(idx_t jblk) const {
        const idx_t jpbegin = jblk * nproma;

        if constexpr(Blocked::rank()==4) {
            const auto* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk = &blocked(jblk, 0, 0, 0);
            auto* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk = &nonblocked(jpbegin, 0, 0);
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_blocked_to_nonblocked_rank4_block<static_nrof()>(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_blocked_to_nonblocked_rank4_block(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
        }
        else if constexpr(Blocked::rank()==3) {
            const auto* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk = &blocked(jblk, 0, 0);
            auto* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk = &nonblocked(jpbegin, 0);
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_blocked_to_nonblocked_rank3_block<static_nrof()>(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_blocked_to_nonblocked_rank3_block(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
        }
        else if constexpr (Blocked::rank()==2) {
            const auto* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk = &blocked(jblk, 0);
            auto* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk = &nonblocked(jpbegin);
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_blocked_to_nonblocked_rank2_block<static_nrof()>(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_blocked_to_nonblocked_rank2_block(raw_blocked_jblk, raw_nonblocked_jblk, nrof);
            }
        }
        else {
            ATLAS_THROW_EXCEPTION("host_copy_blocked_to_nonblocked_mdspan not implemented");
        }
    }

private:
    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_blocked_to_nonblocked_rank4_block(const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                                value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                                [[maybe_unused]] const idx_t nrof) const {
        const idx_t nlev = nonblocked.extent(1);
        const idx_t nvar = nonblocked.extent(2);
        const idx_t point_stride = nlev * nvar;
        const idx_t var_block_stride = nlev * nproma;
        const idx_t lev_block_stride = nproma;
        const idx_t lev_nonblocked_stride = nvar;

        if (loop_order == RelayoutLoopOrder::nproma_innermost) {
            for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jvar = raw_blocked_jblk + jvar * var_block_stride;
                value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jvar = raw_nonblocked_jblk + jvar;
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jvar_jlev = raw_blocked_jblk_jvar + jlev * lev_block_stride;
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jvar_jlev = raw_nonblocked_jblk_jvar + jlev * lev_nonblocked_stride;
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
                    const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * point_stride;
                    for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                        const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof_jvar = raw_blocked_jblk_jrof + jvar * var_block_stride;
                        value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof_jvar = raw_nonblocked_jblk_jrof + jvar;
                        ATLAS_RELAYOUT_SIMD
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            raw_nonblocked_jblk_jrof_jvar[jlev * lev_nonblocked_stride] = raw_blocked_jblk_jrof_jvar[jlev * lev_block_stride];
                        }
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * point_stride;
                    for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                        const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof_jvar = raw_blocked_jblk_jrof + jvar * var_block_stride;
                        value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof_jvar = raw_nonblocked_jblk_jrof + jvar;
                        ATLAS_RELAYOUT_SIMD
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            raw_nonblocked_jblk_jrof_jvar[jlev * lev_nonblocked_stride] = raw_blocked_jblk_jrof_jvar[jlev * lev_block_stride];
                        }
                    }
                }
            }
        }
    }

    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_blocked_to_nonblocked_rank3_block(const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                                value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                                [[maybe_unused]] const idx_t nrof) const {
        const idx_t nlev = nonblocked.extent(1);

        if (loop_order == RelayoutLoopOrder::nproma_innermost) {
            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jlev = raw_blocked_jblk + jlev * nproma;
                value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jlev = raw_nonblocked_jblk + jlev;
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
                    const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * nlev;
                    ATLAS_RELAYOUT_SIMD
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        raw_nonblocked_jblk_jrof[jlev] = raw_blocked_jblk_jrof[jlev * nproma];
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk_jrof = raw_blocked_jblk + jrof;
                    value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk_jrof = raw_nonblocked_jblk + jrof * nlev;
                    ATLAS_RELAYOUT_SIMD
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        raw_nonblocked_jblk_jrof[jlev] = raw_blocked_jblk_jrof[jlev * nproma];
                    }
                }
            }
        }
    }

    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_blocked_to_nonblocked_rank2_block(const value_type* ATLAS_RELAYOUT_RESTRICT raw_blocked_jblk,
                                                value_type* ATLAS_RELAYOUT_RESTRICT raw_nonblocked_jblk,
                                                [[maybe_unused]] const idx_t nrof) const {
        if (use_memcpy) {
            const std::size_t count = static_cast<std::size_t>((nrof_static == 0) ? nrof : nrof_static);
            std::memcpy(raw_nonblocked_jblk, raw_blocked_jblk, count * sizeof(value_type));
        }
        else if constexpr (nrof_static == 0) {
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                raw_nonblocked_jblk[jrof] = raw_blocked_jblk[jrof];
            }
        }
        else {
            for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                raw_nonblocked_jblk[jrof] = raw_blocked_jblk[jrof];
            }
        }
    }
};

template <size_t nproma_extent, class Blocked, class Nonblocked>
auto make_copy_blocked_to_nonblocked_contiguous_raw_pointers(const Blocked blocked, Nonblocked nonblocked) {
    return CopyBlockedToNonblockedContiguousRawPointers<nproma_extent, Blocked, Nonblocked>{blocked, nonblocked};
}

template <size_t nproma_extent, class Blocked, class Nonblocked>
/// @brief Launch the raw-pointer blocked-to-nonblocked kernel for all blocks.
///
/// Preconditions:
/// - `nonblocked.rank() == blocked.rank() - 1`
/// - `blocked` is block-contiguous
/// - `nonblocked` behaves like layout-right contiguous storage
ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
void host_copy_blocked_to_nonblocked_contiguous_raw_pointers_nproma(const Blocked blocked, Nonblocked nonblocked) {
    ATLAS_ASSERT(is_block_contiguous(blocked));
    ATLAS_ASSERT(array::introspection::can_use_layout_right(nonblocked));

    const idx_t nblks = blocked.extent(0);
    const idx_t nproma = last_extent(blocked);
    if constexpr (nproma_extent != dynamic_extent) {
        ATLAS_ASSERT(nproma_extent == nproma);
    }

    static_assert(nonblocked.rank() == blocked.rank()-1);

    auto copy_blocked_to_nonblocked_block = make_copy_blocked_to_nonblocked_contiguous_raw_pointers<nproma_extent>(blocked, nonblocked);
    atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
        copy_blocked_to_nonblocked_block(jblk);
    }
}

#endif

/// @brief Trait family describing one blocked subspan used by the mdspan implementation kernels.
///
/// Specializations adapt the blocked Atlas or mdspan-like view into a layout-right mdspan
/// representing one logical block. `block_alignment` selects aligned or unaligned accessor
/// policies, and `nproma` may be static or dynamic.
template <typename Blocked, size_t nproma, BlockAlignment block_alignment, idx_t Rank = Blocked::rank()>
struct BlockedSubspanTraits;

template <typename Blocked, size_t nproma, BlockAlignment block_alignment>
struct BlockedSubspanTraits<Blocked, nproma, block_alignment, 4> {
    static constexpr bool is_aligned = (block_alignment == BlockAlignment::aligned);
    using value_t = array::introspection::element_t<Blocked>;
    static constexpr std::size_t alignment = (not is_aligned) ? alignof(value_t) : ::atlas::alignment;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t, alignment>;
    using extents_t = extents<idx_t, dynamic_extent, dynamic_extent, nproma>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};

template <typename Blocked, size_t nproma, BlockAlignment block_alignment>
struct BlockedSubspanTraits<Blocked, nproma, block_alignment, 3> {
    static constexpr bool is_aligned = (block_alignment == BlockAlignment::aligned);
    using value_t = array::introspection::element_t<Blocked>;
    static constexpr std::size_t alignment = (not is_aligned) ? alignof(value_t) : ::atlas::alignment;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t, alignment>;
    using extents_t = extents<idx_t, dynamic_extent, nproma>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};

template <typename Blocked, size_t nproma, BlockAlignment block_alignment>
struct BlockedSubspanTraits<Blocked, nproma, block_alignment, 2> {
    static constexpr bool is_aligned = (block_alignment == BlockAlignment::aligned);
    using value_t = array::introspection::element_t<Blocked>;
    static constexpr std::size_t alignment = (nproma == dynamic_extent || not is_aligned) ? alignof(value_t) : blocked_subspan_alignment_v<nproma, value_t>;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t, alignment>;
    using extents_t = extents<idx_t, nproma>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};

template <typename Nonblocked, size_t nproma, idx_t Rank = Nonblocked::rank()>
/// @brief Trait family describing one nonblocked chunk used by the mdspan implementation kernels.
///
/// Each specialization exposes an mdspan type covering the logical point interval that maps
/// to one blocked chunk. These traits let the mdspan implementation kernels operate on regular
/// layout-right mdspan slices irrespective of the original view type.
struct NonblockedSubspanTraits;

template <typename Nonblocked, size_t nproma>
struct NonblockedSubspanTraits<Nonblocked, nproma, 3> {
    using value_t = array::introspection::element_t<Nonblocked>;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t,alignment>;
    using extents_t = extents<idx_t, nproma, dynamic_extent, dynamic_extent>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};

template <typename Nonblocked, size_t nproma>
struct NonblockedSubspanTraits<Nonblocked, nproma, 2> {
    using value_t = array::introspection::element_t<Nonblocked>;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t,alignment>;
    using extents_t = extents<idx_t, nproma, dynamic_extent>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};

template <typename Nonblocked, size_t nproma>
struct NonblockedSubspanTraits<Nonblocked, nproma, 1> {
    using value_t = array::introspection::element_t<Nonblocked>;
    using layout_t = layout_right;
    using accessor_t = restrict_aligned_accessor<value_t,alignment>;
    using extents_t = extents<idx_t, nproma>;
    using span_t = mdspan<value_t, extents_t, layout_t, accessor_t>;
};


template <typename Blocked, size_t nproma_extent = dynamic_extent, BlockAlignment block_alignment = BlockAlignment::aligned>
using blocked_subspan_t = typename BlockedSubspanTraits<Blocked,nproma_extent,block_alignment>::span_t;

template <typename Nonblocked, size_t nproma = dynamic_extent>
using nonblocked_subspan_t = typename NonblockedSubspanTraits<Nonblocked,nproma>::span_t;

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
/// @brief Return whether every block base in a blocked chunk satisfies the alignment contract.
///
/// The mdspan implementation distinguishes aligned and unaligned block accessors. The aligned
/// accessor applies `assume_aligned` to each per-block base pointer `&blocked(jblk, 0, ...)`, so
/// alignment holds for all blocks iff the global base is aligned AND the spacing between
/// consecutive block bases keeps them aligned. That spacing (in elements, allowing for padding
/// between blocks) is `stride(0)`, the stride of the block dimension -- not `stride(1)`, which for
/// a rank-2 view is the innermost `nproma` stride (== 1) and would spuriously report unaligned.
bool is_block_aligned(const Blocked& blocked) {
    using value_type = array::introspection::element_t<decltype(blocked)>;
    return is_aligned(blocked,alignment) &&
           (static_cast<std::size_t>(blocked.stride(0)) * sizeof(value_type) % alignment == 0);
}

template<typename Blocked>
/// @brief Assert the blocked-view preconditions shared by all host relayout kernels.
void assert_requirements_on_blocked(const Blocked& blocked) {
    ATLAS_ASSERT(is_block_contiguous(blocked));
}
template<typename Nonblocked>
/// @brief Assert the nonblocked-view preconditions shared by all host relayout kernels.
void assert_requirements_on_nonblocked(const Nonblocked& nonblocked) {
    bool nonblocked_is_aligned = is_aligned(nonblocked,alignment);
    ATLAS_ASSERT(nonblocked_is_aligned);
    ATLAS_ASSERT(array::introspection::can_use_layout_right(nonblocked));
}

template <size_t nproma_extent, BlockAlignment block_alignment, class Blocked, class Nonblocked>
/// @brief Per-block mdspan implementation kernel for blocked-to-nonblocked host relayout.
///
/// This worker materializes one blocked subspan and one nonblocked chunk as layout-right mdspans,
/// then performs the copy with either aligned or unaligned accessor policies. It is selected when
/// the raw-pointer path is not used or when the mdspan implementation is chosen.
struct CopyBlockedToNonblockedMdspan {
    using blocked_subspan_type = blocked_subspan_t<Blocked, nproma_extent, block_alignment>;
    using nonblocked_subspan_type = nonblocked_subspan_t<Nonblocked, nproma_extent>;
    using block_extents_type = typename blocked_subspan_type::extents_type;
    using nonblocked_extents_type = typename nonblocked_subspan_type::extents_type;

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
    bool use_memcpy;
    block_extents_type block_extents;
    nonblocked_extents_type nonblocked_extents;

    CopyBlockedToNonblockedMdspan(Blocked blocked, Nonblocked nonblocked):
        blocked(blocked),
        nonblocked(nonblocked),
        np(nonblocked.extent(0)),
        nblks(blocked.extent(0)),
        nproma(last_extent(blocked)),
        loop_order(relayout_loop_order()),
        use_memcpy(relayout_blocked_nonblocked_use_memcpy()),
        block_extents(make_blocked_subspan_extents<nproma_extent>(blocked)),
        nonblocked_extents(make_nonblocked_subspan_extents<nproma_extent>(nonblocked)) {}

    /// @brief Copy one blocked chunk into the corresponding nonblocked chunk.
    /// @param jblk Block index in the blocked input view.
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
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
            const auto* blocked_jblk_data_handle = &blocked(jblk, 0);
            auto* nonblocked_jblk_data_handle = &nonblocked(jpbegin);
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_blocked_to_nonblocked_rank2_block<static_nrof()>(blocked_jblk_data_handle, nonblocked_jblk_data_handle, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_blocked_to_nonblocked_rank2_block(blocked_jblk_data_handle, nonblocked_jblk_data_handle, nrof);
            }
        }
        else {
            ATLAS_THROW_EXCEPTION("host_copy_blocked_to_nonblocked_mdspan not implemented");
        }
    }

private:

    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_blocked_to_nonblocked_rank4_block(blocked_subspan_type& block,
                                                const nonblocked_subspan_type& nonblocked_chunk,
                                                [[maybe_unused]] const idx_t nrof) const {
        const idx_t nvar = block.extent(0);
        const idx_t nlev = block.extent(1);
        if constexpr (nrof_static == 0) {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        ATLAS_RELAYOUT_MDSPAN_SIMD
                        for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                            nonblocked_chunk(jrof,jlev,jvar) = block(jvar,jlev,jrof);
                        }
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        ATLAS_RELAYOUT_MDSPAN_SIMD
                        for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                            nonblocked_chunk(jrof,jlev,jvar) = block(jvar,jlev,jrof);
                        }
                    }
                }
            }
        }
        else {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        ATLAS_RELAYOUT_MDSPAN_SIMD
                        for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                            nonblocked_chunk(jrof,jlev,jvar) = block(jvar,jlev,jrof);
                        }
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        ATLAS_RELAYOUT_MDSPAN_SIMD
                        for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                            nonblocked_chunk(jrof,jlev,jvar) = block(jvar,jlev,jrof);
                        }
                    }
                }
            }
        }
    }

    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_blocked_to_nonblocked_rank3_block(blocked_subspan_type& block,
                                                const nonblocked_subspan_type& nonblocked_chunk,
                                                [[maybe_unused]] const idx_t nrof) const {
        const idx_t nlev = block.extent(0);
        if constexpr (nrof_static == 0) {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    ATLAS_RELAYOUT_MDSPAN_SIMD
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        nonblocked_chunk(jrof,jlev) = block(jlev,jrof);
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    ATLAS_RELAYOUT_MDSPAN_SIMD
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        nonblocked_chunk(jrof,jlev) = block(jlev,jrof);
                    }
                }
            }
        }
        else {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    ATLAS_RELAYOUT_MDSPAN_SIMD
                    for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                        nonblocked_chunk(jrof,jlev) = block(jlev,jrof);
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    ATLAS_RELAYOUT_MDSPAN_SIMD
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        nonblocked_chunk(jrof,jlev) = block(jlev,jrof);
                    }
                }
            }
        }
    }

    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_blocked_to_nonblocked_rank2_block(typename blocked_subspan_type::data_handle_type blocked_jblk_data_handle,
                                                typename nonblocked_subspan_type::data_handle_type nonblocked_jblk_data_handle,
                                                [[maybe_unused]] const idx_t nrof) const {
        if (use_memcpy) {
            const std::size_t count = static_cast<std::size_t>((nrof_static == 0) ? nrof : nrof_static);
            std::memcpy(nonblocked_jblk_data_handle, blocked_jblk_data_handle,
                        count * sizeof(typename nonblocked_subspan_type::value_type));
        }
        else {
            blocked_subspan_type block_jblk{blocked_jblk_data_handle, block_extents};
            nonblocked_subspan_type nonblocked_jblk{nonblocked_jblk_data_handle, nonblocked_extents};
            if constexpr (nrof_static == 0) {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    nonblocked_jblk(jrof) = block_jblk(jrof);
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    nonblocked_jblk(jrof) = block_jblk(jrof);
                }
            }
        }
    }
};

template <size_t nproma_extent, BlockAlignment block_alignment, class Blocked, class Nonblocked>
auto make_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    return CopyBlockedToNonblockedMdspan<nproma_extent, block_alignment, Blocked, Nonblocked>{blocked, nonblocked};
}

template <size_t nproma_extent, BlockAlignment block_alignment, class Nonblocked, class Blocked>
/// @brief Per-block mdspan implementation kernel for nonblocked-to-blocked host relayout.
///
/// This is the inverse of `CopyBlockedToNonblockedMdspan`. It presents the logical source
/// chunk and destination block as mdspan slices and copies one block at a time while preserving
/// the selected loop order and aligned/unaligned accessor policy.
struct CopyNonblockedToBlockedMdspan {
    using blocked_subspan_type = blocked_subspan_t<Blocked, nproma_extent, block_alignment>;
    using nonblocked_subspan_type = nonblocked_subspan_t<Nonblocked, nproma_extent>;
    using block_extents_type = typename blocked_subspan_type::extents_type;
    using nonblocked_extents_type = typename nonblocked_subspan_type::extents_type;

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
    bool use_memcpy;
    block_extents_type block_extents;
    nonblocked_extents_type nonblocked_extents;

    CopyNonblockedToBlockedMdspan(Nonblocked nonblocked, Blocked blocked):
        nonblocked(nonblocked),
        blocked(blocked),
        np(nonblocked.extent(0)),
        nblks(blocked.extent(0)),
        nproma(last_extent(blocked)),
        loop_order(relayout_loop_order()),
        use_memcpy(relayout_blocked_nonblocked_use_memcpy()),
        block_extents(make_blocked_subspan_extents<nproma_extent>(blocked)),
        nonblocked_extents(make_nonblocked_subspan_extents<nproma_extent>(nonblocked)) {}

    /// @brief Copy one nonblocked chunk into the corresponding blocked chunk.
    /// @param jblk Block index in the blocked output view.
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
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
            const auto* nonblocked_jblk_data_handle = &nonblocked(jpbegin);
            auto* block_jblk_data_handle = &blocked(jblk, 0);
            if (jblk < nblks-1) {
                const idx_t nrof = nproma;
                copy_nonblocked_to_blocked_rank2_block<static_nrof()>(nonblocked_jblk_data_handle, block_jblk_data_handle, nrof);
            }
            else {
                const idx_t nrof = std::min(np - jpbegin, nproma);
                copy_nonblocked_to_blocked_rank2_block(nonblocked_jblk_data_handle, block_jblk_data_handle, nrof);
            }
        }
        else {
            ATLAS_THROW_EXCEPTION("host_copy_nonblocked_to_blocked_mdspan not implemented");
        }
    }

private:

    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_nonblocked_to_blocked_rank4_block(const nonblocked_subspan_type& nonblocked_chunk, blocked_subspan_type& block,
                                                [[maybe_unused]] const idx_t nrof) const {
        const idx_t nvar = block.extent(0);
        const idx_t nlev = block.extent(1);
        if constexpr (nrof_static == 0) {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        ATLAS_RELAYOUT_MDSPAN_SIMD
                        for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                            block(jvar,jlev,jrof) = nonblocked_chunk(jrof,jlev,jvar);
                        }
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        ATLAS_RELAYOUT_MDSPAN_SIMD
                        for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                            block(jvar,jlev,jrof) = nonblocked_chunk(jrof,jlev,jvar);
                        }
                    }
                }
            }
        }
        else {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        ATLAS_RELAYOUT_MDSPAN_SIMD
                        for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                            block(jvar,jlev,jrof) = nonblocked_chunk(jrof,jlev,jvar);
                        }
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        ATLAS_RELAYOUT_MDSPAN_SIMD
                        for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                            block(jvar,jlev,jrof) = nonblocked_chunk(jrof,jlev,jvar);
                        }
                    }
                }
            }
        }
    }

    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_nonblocked_to_blocked_rank3_block(const nonblocked_subspan_type& nonblocked_chunk, blocked_subspan_type& block,
                                                [[maybe_unused]] const idx_t nrof) const {
        const idx_t nlev = block.extent(0);
        if constexpr (nrof_static == 0) {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    ATLAS_RELAYOUT_MDSPAN_SIMD
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        block(jlev,jrof) = nonblocked_chunk(jrof,jlev);
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    ATLAS_RELAYOUT_MDSPAN_SIMD
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        block(jlev,jrof) = nonblocked_chunk(jrof,jlev);
                    }
                }
            }
        }
        else {
            if (loop_order == RelayoutLoopOrder::nproma_innermost) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    ATLAS_RELAYOUT_MDSPAN_SIMD
                    for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                        block(jlev,jrof) = nonblocked_chunk(jrof,jlev);
                    }
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    ATLAS_RELAYOUT_MDSPAN_SIMD
                    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                        block(jlev,jrof) = nonblocked_chunk(jrof,jlev);
                    }
                }
            }
        }
    }

    template <idx_t nrof_static = 0>
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_nonblocked_to_blocked_rank2_block(typename nonblocked_subspan_type::data_handle_type nonblocked_jblk_data_handle,
                                                typename blocked_subspan_type::data_handle_type block_jblk_data_handle,
                                                [[maybe_unused]] const idx_t nrof) const {
        if (use_memcpy) {
            const std::size_t count = static_cast<std::size_t>((nrof_static == 0) ? nrof : nrof_static);
            std::memcpy(block_jblk_data_handle, nonblocked_jblk_data_handle,
                        count * sizeof(typename blocked_subspan_type::value_type));
        }
        else {
            blocked_subspan_type block_jblk{block_jblk_data_handle, block_extents};
            nonblocked_subspan_type nonblocked_jblk{nonblocked_jblk_data_handle, nonblocked_extents};
            if constexpr (nrof_static == 0) {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    block_jblk(jrof) = nonblocked_jblk(jrof);
                }
            }
            else {
                for (idx_t jrof = 0; jrof < nrof_static; ++jrof) {
                    block_jblk(jrof) = nonblocked_jblk(jrof);
                }
            }
        }
    }

};

template <size_t nproma_extent, BlockAlignment block_alignment, class Nonblocked, class Blocked>
auto make_copy_nonblocked_to_blocked_block(const Nonblocked nonblocked, Blocked blocked) {
    return CopyNonblockedToBlockedMdspan<nproma_extent, block_alignment, Nonblocked, Blocked>{nonblocked, blocked};
}

template <size_t nproma_extent, class Blocked, class Nonblocked>
/// @brief Run the mdspan implementation blocked-to-nonblocked kernel for a fixed or dynamic `nproma`.
///
/// This wrapper validates the common preconditions, chooses aligned or unaligned block access,
/// and then launches one per-block worker over the blocked extent.
ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
void host_copy_blocked_to_nonblocked_nproma(const Blocked blocked, Nonblocked nonblocked) {
    assert_requirements_on_blocked(blocked);
    assert_requirements_on_nonblocked(nonblocked);

    // Add a harmless side-effect to force a unique assembly signature
    volatile int linker_poison = 42; 
    (void)linker_poison; 

    const idx_t nblks = blocked.extent(0);
    if (is_block_aligned(blocked)) {
        auto copy_blocked_to_nonblocked_block = make_copy_blocked_to_nonblocked_mdspan<nproma_extent, BlockAlignment::aligned>(blocked, nonblocked);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            copy_blocked_to_nonblocked_block(jblk);
        }
    }
    else {
        Log::debug() << "host_copy_blocked_to_nonblocked_nproma: A block is not aligned, using the unaligned mdspan implementation.";
        Log::debug() << "\nBlocked: " << std::vector<idx_t>(blocked.shape(), blocked.shape() + blocked.rank()) << std::endl;

        auto copy_blocked_to_nonblocked_block = make_copy_blocked_to_nonblocked_mdspan<nproma_extent, BlockAlignment::unaligned>(blocked, nonblocked);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            copy_blocked_to_nonblocked_block(jblk);
        }
    }
}

template <size_t nproma_extent, class Nonblocked, class Blocked>
/// @brief Run the mdspan implementation nonblocked-to-blocked kernel for a fixed or dynamic `nproma`.
///
/// This wrapper validates the common preconditions, chooses aligned or unaligned block access,
/// and then launches one per-block worker over the blocked extent.
ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
void host_copy_nonblocked_to_blocked_nproma(const Nonblocked nonblocked, Blocked blocked) {
    assert_requirements_on_nonblocked(nonblocked);
    assert_requirements_on_blocked(blocked);

    const idx_t nblks = blocked.extent(0);
    if (is_block_aligned(blocked)) {
        auto copy_nonblocked_to_blocked_block = make_copy_nonblocked_to_blocked_block<nproma_extent, BlockAlignment::aligned>(nonblocked, blocked);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            copy_nonblocked_to_blocked_block(jblk);
        }
    }
    else {
        Log::debug() << "host_copy_nonblocked_to_blocked_nproma: A block is not aligned, using the unaligned mdspan implementation.";
        Log::debug() << "\nBlocked: " << std::vector<idx_t>(blocked.shape(), blocked.shape() + blocked.rank()) << std::endl;

        auto copy_nonblocked_to_blocked_block = make_copy_nonblocked_to_blocked_block<nproma_extent, BlockAlignment::unaligned>(nonblocked, blocked);
        atlas_omp_parallel_for(idx_t jblk = 0; jblk < nblks; ++jblk) {
            copy_nonblocked_to_blocked_block(jblk);
        }
    }
}

}  // namespace

template <class Nonblocked, class Blocked>
/// @brief Dispatch nonblocked-to-blocked host relayout to the best available kernel family.
///
/// Selection order:
/// - raw-pointer contiguous kernel when implementation mode is `raw_pointers` and both layouts
///   allow direct pointer access;
/// - otherwise the mdspan/subspan implementation;
/// - within each family, prefer compile-time `nproma` specialization when available, otherwise
///   switch over common runtime `nproma` values before using `dynamic_extent`.
ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
void host_copy_nonblocked_to_blocked_impl(const Nonblocked nonblocked, Blocked blocked) {
    static_assert(nonblocked.rank() == blocked.rank()-1);
    const idx_t nproma = blocked.extent(blocked.rank() - 1);

    // If Blocked is an mdspan whose last (nproma) dimension is a static extent, then nproma is
    // known at compile time.  In that case we can call the statically-sized kernels directly and
    // skip the runtime switch dispatch below, avoiding the unused template instantiations.
    constexpr std::size_t nproma_extent = nproma_extent_v<Blocked>;

    ATLAS_ASSERT(array::introspection::can_use_layout_right(nonblocked));
    ATLAS_ASSERT(is_block_contiguous(blocked));

#if DISABLE_RAW_POINTERS == 0
    if (relayout_implementation() == RelayoutImplementation::raw_pointers) {
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

    // Use the mdspan implementation when the raw-pointer implementation is not selected.
    // This path still requires a contiguous nonblocked view and a block-contiguous blocked view,
    // enforced by `assert_requirements_on_nonblocked()` and `assert_requirements_on_blocked()`.
    // Apply static dispatch there as well for better performance.
    if (relayout_nproma_dispatch() != RelayoutNpromaDispatch::static_dispatch) {
        // Runtime dispatch for the mdspan implementation
        host_copy_nonblocked_to_blocked_nproma<dynamic_extent>(nonblocked, blocked);
    }

    // Optimized mdspan implementation paths with static nproma dispatch.
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
 * @brief Dispatch blocked-to-nonblocked host relayout to the best available kernel family.
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
 *
 * Dispatch order mirrors `host_copy_nonblocked_to_blocked_impl`: choose the raw-pointer fast
 * path when possible, otherwise use the mdspan/subspan implementation, and in both cases prefer
 * statically-sized `nproma` instantiations over runtime-generic variants.
 */
void host_copy_blocked_to_nonblocked_impl(const Blocked blocked, Nonblocked nonblocked) {
    static_assert(nonblocked.rank() == blocked.rank()-1);

    // If Blocked is an mdspan whose last (nproma) dimension is a static extent, then nproma is
    // known at compile time.  In that case we can call the statically-sized kernels directly and
    // skip the runtime switch dispatch below, avoiding the unused template instantiations.
    constexpr std::size_t nproma_extent = nproma_extent_v<Blocked>;

    ATLAS_ASSERT(array::introspection::can_use_layout_right(nonblocked));
    ATLAS_ASSERT(is_block_contiguous(blocked));

#if DISABLE_RAW_POINTERS == 0
    if (relayout_implementation() == RelayoutImplementation::raw_pointers) {
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

    // Use the mdspan implementation when the raw-pointer implementation is not selected.
    // This path still requires a contiguous nonblocked view and a block-contiguous blocked view,
    // enforced by `assert_requirements_on_nonblocked()` and `assert_requirements_on_blocked()`.
    // Apply static dispatch there as well for better performance.
    if (relayout_nproma_dispatch() != RelayoutNpromaDispatch::static_dispatch) {
        // Runtime dispatch for the mdspan implementation
        return host_copy_blocked_to_nonblocked_nproma<dynamic_extent>(blocked, nonblocked);
    }

    // Optimized mdspan implementation paths with static nproma dispatch.
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
/// @brief Per-block raw-pointer kernel for blocked-to-blocked host relayout.
///
/// This worker supports different input and output `nproma` values by walking the logical point
/// range and copying maximal contiguous chunks between the source block and destination block.
/// Rank-specific helpers preserve the inner non-horizontal layout while remapping the horizontal
/// blocked dimension.
struct CopyBlockedToBlockedBlockRawPointers {
    using value_t = std::decay_t<typename BlockedOut::value_type>;

    const BlockedIn blocked_in;
    mutable BlockedOut blocked_out;
    idx_t nblks_in;
    idx_t nproma_in;
    idx_t nblks_out;
    idx_t nproma_out;
    idx_t total_points;
    bool use_memcpy;

    CopyBlockedToBlockedBlockRawPointers(const BlockedIn blocked_in, BlockedOut blocked_out):
        blocked_in(blocked_in),
        blocked_out(blocked_out),
        nblks_in(blocked_in.extent(0)),
        nproma_in(blocked_in.extent(blocked_in.rank()-1)),
        nblks_out(blocked_out.extent(0)),
        nproma_out(blocked_out.extent(blocked_out.rank()-1)),
        total_points(std::min(nblks_in * nproma_in, nblks_out * nproma_out)),
        use_memcpy(relayout_blocked_to_blocked_use_memcpy()) {
        static_assert(std::is_same_v<std::decay_t<typename BlockedIn::value_type>, std::decay_t<typename BlockedOut::value_type>>, "Data types of input and output views must match for blocked-to-blocked copy");
        static_assert(BlockedIn::rank() == BlockedOut::rank());

        if constexpr (BlockedIn::rank() == 4) {
            ATLAS_ASSERT(blocked_in.extent(1) == blocked_out.extent(1));
            ATLAS_ASSERT(blocked_in.extent(2) == blocked_out.extent(2));
        }
        else if constexpr (BlockedIn::rank() == 3) {
            ATLAS_ASSERT(blocked_in.extent(1) == blocked_out.extent(1));
        }
    }

    /// @brief Copy one destination block from the corresponding logical source range.
    /// @param jblk_out Block index in the blocked output view.
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void operator()(idx_t jblk_out) const {
        const idx_t jpbegin = jblk_out * nproma_out;
        if (jpbegin >= total_points) {
            return;
        }

        if constexpr (BlockedIn::rank()==4) {
            value_t* ATLAS_RELAYOUT_RESTRICT raw_blocked_out_jblk = &blocked_out(jblk_out, 0, 0, 0);
            copy_rank4_block(raw_blocked_out_jblk, jblk_out, jpbegin);
        }
        else if constexpr (BlockedIn::rank()==3) {
            value_t* ATLAS_RELAYOUT_RESTRICT raw_blocked_out_jblk = &blocked_out(jblk_out, 0, 0);
            copy_rank3_block(raw_blocked_out_jblk, jblk_out, jpbegin);
        }
        else if constexpr (BlockedIn::rank()==2) {
            value_t* ATLAS_RELAYOUT_RESTRICT raw_blocked_out_jblk = &blocked_out(jblk_out, 0);
            copy_rank2_block(raw_blocked_out_jblk, jblk_out, jpbegin);
        }
        else {
            ATLAS_THROW_EXCEPTION("transposition not implemented for rank " << blocked_in.rank());
        }
    }

private:
    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_rank4_block(value_t* ATLAS_RELAYOUT_RESTRICT raw_blocked_out_jblk,
                          const idx_t jblk_out,
                          const idx_t jpbegin) const {
        const idx_t nlev = blocked_in.extent(1);
        const idx_t nvar = blocked_in.extent(2);
        const idx_t out_lev_stride = nvar * nproma_out;
        const idx_t out_var_stride = nproma_out;
        const idx_t in_lev_stride = nvar * nproma_in;
        const idx_t in_var_stride = nproma_in;
        const idx_t jpend = std::min(total_points, jpbegin + nproma_out);
        idx_t jp = jpbegin;

        while (jp < jpend) {
            const idx_t jblk_in  = jp / nproma_in;
            const idx_t jrof_in  = jp - jblk_in * nproma_in;
            const idx_t jrof_out = jp - jblk_out * nproma_out;
            const idx_t chunk = std::min(jpend - jp, nproma_in - jrof_in);
            const value_t* ATLAS_RELAYOUT_RESTRICT raw_blocked_in_jblk = &blocked_in(jblk_in, 0, 0, 0);

            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                    const idx_t index_out_base = jlev * out_lev_stride + jvar * out_var_stride + jrof_out;
                    const idx_t index_in_base = jlev * in_lev_stride + jvar * in_var_stride + jrof_in;
                    copy_chunk(raw_blocked_out_jblk + index_out_base, raw_blocked_in_jblk + index_in_base, chunk);
                }
            }
            jp += chunk;
        }
    }

    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_rank3_block(value_t* ATLAS_RELAYOUT_RESTRICT raw_blocked_out_jblk,
                          const idx_t jblk_out,
                          const idx_t jpbegin) const {
        const idx_t nlev = blocked_in.extent(1);
        const idx_t out_lev_stride = nproma_out;
        const idx_t in_lev_stride = nproma_in;
        const idx_t jpend = std::min(total_points, jpbegin + nproma_out);
        idx_t jp = jpbegin;

        while (jp < jpend) {
            const idx_t jblk_in  = jp / nproma_in;
            const idx_t jrof_in  = jp - jblk_in * nproma_in;
            const idx_t jrof_out = jp - jblk_out * nproma_out;
            const idx_t chunk = std::min(jpend - jp, nproma_in - jrof_in);
            const value_t* ATLAS_RELAYOUT_RESTRICT raw_blocked_in_jblk = &blocked_in(jblk_in, 0, 0);

            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                const idx_t index_out_base = jlev * out_lev_stride + jrof_out;
                const idx_t index_in_base = jlev * in_lev_stride + jrof_in;
                copy_chunk(raw_blocked_out_jblk + index_out_base, raw_blocked_in_jblk + index_in_base, chunk);
            }
            jp += chunk;
        }
    }

    ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
    void copy_rank2_block(value_t* ATLAS_RELAYOUT_RESTRICT raw_blocked_out_jblk,
                          const idx_t jblk_out,
                          const idx_t jpbegin) const {
        const idx_t jpend = std::min(total_points, jpbegin + nproma_out);
        idx_t jp = jpbegin;

        while (jp < jpend) {
            const idx_t jblk_in  = jp / nproma_in;
            const idx_t jrof_in  = jp - jblk_in * nproma_in;
            const idx_t jrof_out = jp - jblk_out * nproma_out;
            const idx_t chunk = std::min(jpend - jp, nproma_in - jrof_in);
            const value_t* ATLAS_RELAYOUT_RESTRICT raw_blocked_in_jblk = &blocked_in(jblk_in, 0);
            copy_chunk(raw_blocked_out_jblk + jrof_out, raw_blocked_in_jblk + jrof_in, chunk);
            jp += chunk;
        }
    }

    void copy_chunk(value_t* ATLAS_RELAYOUT_RESTRICT out,
                    const value_t* ATLAS_RELAYOUT_RESTRICT in,
                    const idx_t chunk) const {
        if (use_memcpy) {
            std::memcpy(out, in, static_cast<std::size_t>(chunk) * sizeof(value_t));
        }
        else {
            ATLAS_RELAYOUT_SIMD
            for (idx_t j = 0; j < chunk; ++j) {
                out[j] = in[j];
            }
        }
    }
};

template <class BlockedIn, class BlockedOut>
auto make_copy_blocked_to_blocked(const BlockedIn blocked_in, BlockedOut blocked_out) {
    return CopyBlockedToBlockedBlockRawPointers<BlockedIn, BlockedOut>{blocked_in, blocked_out};
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
 *
 * When enabled and safe, the implementation first attempts a single whole-view `memcpy` for the
 * equal-`nproma`, fully contiguous case. Otherwise it falls back to the per-block raw-pointer
 * kernel, which preserves logical ordering while repacking points into the destination block size.
 */
ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
void host_copy_blocked_to_blocked_impl(const BlockedIn blocked_in, BlockedOut blocked_out) {
    static_assert(std::is_same_v<std::decay_t<typename BlockedIn::value_type>, std::decay_t<typename BlockedOut::value_type>>, "Data types of input and output views must match for blocked-to-blocked copy");
    using value_type = std::decay_t<typename BlockedOut::value_type>;
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

    // At the moment we implement only contiguous block copies for optimizations
    ATLAS_ASSERT(array::introspection::can_use_layout_right(blocked_in));
    ATLAS_ASSERT(array::introspection::can_use_layout_right(blocked_out));

    const bool use_memcpy = relayout_blocked_to_blocked_use_memcpy();

    if (use_memcpy) {
        if (nproma_in == nproma_out && blocked_in.size() == blocked_out.size()) {
            const value_type* raw_in = array::introspection::data_handle(blocked_in);
            value_type* raw_out = array::introspection::data_handle(blocked_out);
            std::memcpy(raw_out, raw_in, blocked_out.size() * sizeof(value_type));
            return;
        }
    }

    auto copy_blocked_to_blocked_block = make_copy_blocked_to_blocked(blocked_in, blocked_out);
    atlas_omp_parallel_for(idx_t jblk_out = 0; jblk_out < nblks_out; ++jblk_out) {
        copy_blocked_to_blocked_block(jblk_out);
    }
}

template <class Nonblocked, class Blocked>
ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
void host_copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked) {
    return host_copy_nonblocked_to_blocked_impl(nonblocked, blocked);
}

template <class Blocked, class Nonblocked>
ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
void host_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    return host_copy_blocked_to_nonblocked_impl(blocked, nonblocked);
}

template <class BlockedIn, class BlockedOut>
ATLAS_RELAYOUT_NOINLINE_IF_PROFILING
void host_copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out) {
    return host_copy_blocked_to_blocked_impl(blocked_in, blocked_out);
}

#undef ATLAS_RELAYOUT_RESTRICT

}  // namespace atlas


#define ATLAS_RELAYOUT_EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(TYPE, BLOCKED_RANK) \
    template void atlas::host_copy_blocked_to_nonblocked_mdspan<atlas::array::ArrayView<const TYPE,BLOCKED_RANK>,atlas::array::ArrayView<TYPE,BLOCKED_RANK-1>>(atlas::array::ArrayView<const TYPE,BLOCKED_RANK>, atlas::array::ArrayView<TYPE,BLOCKED_RANK-1>); \
    template void atlas::host_copy_nonblocked_to_blocked_mdspan<atlas::array::ArrayView<const TYPE,BLOCKED_RANK-1>,atlas::array::ArrayView<TYPE,BLOCKED_RANK>>(atlas::array::ArrayView<const TYPE,BLOCKED_RANK-1>, atlas::array::ArrayView<TYPE,BLOCKED_RANK>); \
    template void atlas::host_copy_blocked_to_blocked_mdspan<atlas::array::ArrayView<const TYPE,BLOCKED_RANK>,atlas::array::ArrayView<TYPE,BLOCKED_RANK>>(atlas::array::ArrayView<const TYPE,BLOCKED_RANK>, atlas::array::ArrayView<TYPE,BLOCKED_RANK>);
