/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/runtime/Exception.h"

#include <cstdlib>
#include <cstring>
#include <type_traits>

/**
 * @file relayout_on_device.tcc
 * @brief Device implementations for copying Atlas data between blocked and nonblocked layouts.
 *
 * Device relayout expects device-accessible Atlas views (`atlas::View`/`atlas::ArrayView`) or
 * mdspan-like views with the rank, datatype, and shape requirements documented in
 * `atlas/util/relayout.h`.
 */

#include "hic/hic.h"
#include "atlas/array.h"
#include "atlas/mdspan.h"
#include "atlas/util/relayout.h"

#include "atlas/runtime/Log.h"

namespace atlas {

ATLAS_DEVICE inline idx_t min(idx_t a, idx_t b) {
    return a < b ? a : b;
}

class KernelTraversal {
public:
    idx_t offset;
    idx_t stride;

    ATLAS_HOST_DEVICE static KernelTraversal current() {
        KernelTraversal traversal;
    #if HIC_COMPILER
        traversal.offset = static_cast<idx_t>(blockIdx.x * blockDim.x + threadIdx.x);
        traversal.stride = static_cast<idx_t>(blockDim.x * gridDim.x);
    #else
        traversal.offset = 0;
        traversal.stride = 1;
    #endif
        return traversal;
    }
};

namespace {
constexpr const char* relayout_implementation_env = "ATLAS_RELAYOUT_IMPLEMENTATION";

enum class RelayoutImplementation {
    arrayview,
    mdspan
};

RelayoutImplementation relayout_implementation() {
    static RelayoutImplementation cached = []() {
         const char* impl = std::getenv(relayout_implementation_env);
         if (impl && std::strcmp(impl, "mdspan") == 0) {
             return RelayoutImplementation::mdspan;
         }
         return RelayoutImplementation::arrayview;
    }();
    return cached;
}

constexpr const char* relayout_loop_order_env = "ATLAS_RELAYOUT_LOOP_ORDER";

// Device relayout loop order, selected via ATLAS_RELAYOUT_LOOP_ORDER (shared with the host
// relayout and the atlas-benchmark-relayout --loop-order option).
//
//  - nproma_outermost (device default): the flat, element-per-thread grid-stride kernels.
//    Consecutive threads map to consecutive points, which keeps the contiguous (nproma)
//    dimension coalesced on the blocked side. This preserves the historical device behaviour.
//  - nproma_innermost: the per-block kernels. Each work item owns one blocked chunk and loops the
//    contiguous nproma dimension as the innermost sequential loop, hoisting the per-element
//    div/mod. This mirrors the host nproma_innermost loop order and is offered so the tradeoff
//    (contiguous per-thread runs vs. cross-thread coalescing) can be measured on real hardware.
//  - nonblocked_coalesced: maps consecutive threads to the innermost (unit-stride) dimension of the
//    nonblocked view instead of to the horizontal point index, so the nonblocked side is coalesced
//    and the blocked side becomes the scattered side (striding by nproma, a small stride, rather
//    than by nlev). For blocked->nonblocked this coalesces the nonblocked writes; for
//    nonblocked->blocked it coalesces the nonblocked reads. Both directions keep the scattered
//    blocked accesses tight (nproma-strided), which can outperform nproma_outermost. The
//    blocked->blocked direction has no nonblocked view and behaves like nproma_outermost.
//  - coalesced_write / coalesced_read: direction-aware meta-orders that select, per direction,
//    whichever concrete kernel above coalesces the write (respectively the read) side. Because the
//    blocked side is a write for nonblocked->blocked but a read for blocked->nonblocked, a fixed
//    kernel cannot be "the coalesced-write kernel" for both directions; these meta-orders resolve
//    to the right kernel for each direction (see resolve_loop_order). Empirically, coalescing the
//    write side is the better default, so coalesced_write resolves to the fastest kernel measured
//    for each direction.
enum class RelayoutLoopOrder {
    nproma_innermost,
    nproma_outermost,
    nonblocked_coalesced,
    coalesced_write,
    coalesced_read
};

RelayoutLoopOrder relayout_loop_order() {
    static RelayoutLoopOrder cached = []() {
         const char* order = std::getenv(relayout_loop_order_env);
         if (order && std::strcmp(order, "nproma_innermost") == 0) {
             return RelayoutLoopOrder::nproma_innermost;
         }
         if (order && std::strcmp(order, "nonblocked_coalesced") == 0) {
             return RelayoutLoopOrder::nonblocked_coalesced;
         }
         if (order && std::strcmp(order, "coalesced_write") == 0) {
             return RelayoutLoopOrder::coalesced_write;
         }
         if (order && std::strcmp(order, "coalesced_read") == 0) {
             return RelayoutLoopOrder::coalesced_read;
         }
         return RelayoutLoopOrder::nproma_outermost;
    }();
    return cached;
}

// The copy direction a kernel implements, used to resolve the direction-aware coalesced_write /
// coalesced_read meta-orders to a concrete kernel choice. The blocked side is the write for
// nonblocked->blocked and the read for blocked->nonblocked, so the same "coalesce the write"
// request maps to different kernels per direction.
enum class CopyDirection {
    nonblocked_to_blocked,
    blocked_to_nonblocked,
    blocked_to_blocked
};

// Resolve coalesced_write / coalesced_read into one of the concrete kernel-selecting loop orders
// (nproma_outermost coalesces the blocked side; nonblocked_coalesced coalesces the nonblocked
// side). Non-meta orders are returned unchanged.
RelayoutLoopOrder resolve_loop_order(RelayoutLoopOrder order, CopyDirection direction) {
    const bool want_write = (order == RelayoutLoopOrder::coalesced_write);
    const bool want_read  = (order == RelayoutLoopOrder::coalesced_read);
    if (!want_write && !want_read) {
        return order;
    }
    switch (direction) {
        case CopyDirection::nonblocked_to_blocked:
            // blocked side = write, nonblocked side = read.
            return want_write ? RelayoutLoopOrder::nproma_outermost      // coalesce blocked write
                              : RelayoutLoopOrder::nonblocked_coalesced; // coalesce nonblocked read
        case CopyDirection::blocked_to_nonblocked:
            // blocked side = read, nonblocked side = write.
            return want_write ? RelayoutLoopOrder::nonblocked_coalesced  // coalesce nonblocked write
                              : RelayoutLoopOrder::nproma_outermost;     // coalesce blocked read
        case CopyDirection::blocked_to_blocked:
            // No nonblocked view; the flat kernel coalesces both blocked read and write.
            return RelayoutLoopOrder::nproma_outermost;
    }
    return RelayoutLoopOrder::nproma_outermost;
}

// Distinguishes the two roles a view can play in a relayout, used as an explicit dispatch tag
// instead of a bare bool. A Blocked view has a contiguous block as its last (nproma) dimension;
// a Nonblocked view is a plain layout_right field.
enum class ViewType {
    Blocked,
    Nonblocked
};

template <typename Layout, typename View>
auto make_relayout_mdspan(View& view) {
    return make_mdspan<Layout, restrict_accessor>(view);
}

template <typename View>
auto make_relayout_mdspan(View& view) {
    return make_relayout_mdspan<layout_stride>(view);
}

template <typename Layout, ViewType VType, typename View, typename Operation>
void dispatch_relayout_view_mdspan(View& view, Operation&& operation) {
    using namespace array::introspection;
    if (is_aligned<64>(view)) {
        if constexpr (VType == ViewType::Blocked) {
            switch (last_extent(view)) {
                case 16: return operation(make_mdspan<extent_with_static_last_dim<16>, Layout, restrict_aligned_accessor_policy<64>>(view));
                case 32: return operation(make_mdspan<extent_with_static_last_dim<32>, Layout, restrict_aligned_accessor_policy<64>>(view));
                case 64: return operation(make_mdspan<extent_with_static_last_dim<64>, Layout, restrict_aligned_accessor_policy<64>>(view));
                default: return operation(make_mdspan<Layout, restrict_aligned_accessor_policy<64>>(view));
            }
        }
        return operation(make_mdspan<Layout, restrict_aligned_accessor_policy<64>>(view));
    }

    if constexpr (VType == ViewType::Blocked) {
        switch (last_extent(view)) {
            case 16: return operation(make_mdspan<extent_with_static_last_dim<16>, Layout, restrict_accessor>(view));
            case 32: return operation(make_mdspan<extent_with_static_last_dim<32>, Layout, restrict_accessor>(view));
            case 64: return operation(make_mdspan<extent_with_static_last_dim<64>, Layout, restrict_accessor>(view));
            default: return operation(make_relayout_mdspan<Layout>(view));
        }
    }
    return operation(make_relayout_mdspan<Layout>(view));
}

template <ViewType Source, ViewType Target, typename SourceView, typename TargetView, typename Operation>
void dispatch_relayout_mdspan(SourceView& source, TargetView& target, Operation&& operation) {
    // Nonblocked fields are always layout_right. Blocked fields always have a contiguous block
    // (the last dimension), but the field as a whole may not be contiguous, so they use
    // layout_stride. The layout for each view is therefore fixed by whether it is blocked.
    using SourceLayout = std::conditional_t<Source == ViewType::Blocked, layout_stride, layout_right>;
    using TargetLayout = std::conditional_t<Target == ViewType::Blocked, layout_stride, layout_right>;
    return dispatch_relayout_view_mdspan<SourceLayout, Source>(source, [&](const auto& source_view) {
        return dispatch_relayout_view_mdspan<TargetLayout, Target>(target, [&](const auto& target_view) {
            return operation(source_view, target_view);
        });
    });
}

#if HIC_COMPILER
// These relayout kernels are simple memory-bound copies with light index arithmetic.
// A 256-thread block is a conservative default on CUDA/HIP GPUs: it is large enough
// to expose memory parallelism, while usually avoiding the register/occupancy tradeoffs
// that can appear with larger 512- or 1024-thread blocks.
constexpr int threads_per_block() {
    return 256;
}

// Grid-stride loops do not need one block per chunk of work. We cap the launch at a
// small multiple of the device SM count so each SM sees several resident blocks while
// the loop carries the remaining work. Four blocks per SM is a pragmatic default for
// these bandwidth-oriented kernels and keeps the heuristic explicit and easy to tune.
constexpr int blocks_per_sm() {
    return 4;
}

int sm_count() {
    int device = 0;
    HIC_CALL(hicGetDevice(&device));

    hicDeviceProp_t properties;
    HIC_CALL(hicGetDeviceProperties(&properties, device));

    return properties.multiProcessorCount;
}

/**
 * @brief Compute the number of blocks to launch for a 1-D grid-stride kernel.
 *
 * The launch uses ceil(work_size / threads_per_block()) blocks, capped to
 * blocks_per_sm() * sm_count(). The cap avoids oversubscribing tiny kernels with
 * very large grids, while grid-stride looping still guarantees full coverage of
 * all work items.
 */
int blocks_per_grid(idx_t work_size) {
    const int max_blocks = blocks_per_sm() * sm_count();
    const auto blocks = static_cast<int>((work_size + threads_per_block() - 1) / threads_per_block());
    return blocks < max_blocks ? blocks : max_blocks;
}
#endif
}

template <class Blocked, class Nonblocked>
ATLAS_GLOBAL void kernel_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    auto npts   = nonblocked.extent(0);
    auto nproma = blocked.extent(blocked.rank()-1);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        idx_t work_size = npts * nlev * nvar;
        auto traversal = KernelTraversal::current();
        for (idx_t offset = traversal.offset; offset < work_size; offset += traversal.stride) {
            idx_t point = offset % npts;
            idx_t entry = offset / npts;
            idx_t jlev  = entry % nlev;
            idx_t jvar  = entry / nlev;
            idx_t jblk  = point / nproma;
            idx_t jrof  = point - jblk * nproma;
            nonblocked(point, jlev, jvar) = blocked(jblk, jvar, jlev, jrof);
        }
    }
    else if constexpr (blocked.rank()==3) {
        idx_t nlev = nonblocked.extent(1);
        idx_t work_size = npts * nlev;
        auto traversal = KernelTraversal::current();
        for (idx_t offset = traversal.offset; offset < work_size; offset += traversal.stride) {
            idx_t point = offset % npts;
            idx_t jlev  = offset / npts;
            idx_t jblk  = point / nproma;
            idx_t jrof  = point - jblk * nproma;
            nonblocked(point, jlev) = blocked(jblk, jlev, jrof);
        }
    }
    else if constexpr (blocked.rank()==2) {
        auto traversal = KernelTraversal::current();
        for (idx_t point = traversal.offset; point < npts; point += traversal.stride) {
            idx_t jblk = point / nproma;
            idx_t jrof = point - jblk * nproma;
            nonblocked(point) = blocked(jblk, jrof);
        }
    }
}

/// @brief Per-block (nproma_innermost) variant of the blocked-to-nonblocked kernel.
///
/// Each grid-stride work item owns one blocked chunk (jblk) for a fixed (jlev, jvar) and loops over
/// the contiguous nproma dimension (jrof) as the innermost sequential loop. This hoists the
/// per-element div/mod out of the inner loop and makes the blocked-side access unit-stride within a
/// thread, mirroring the host nproma_innermost loop order.
template <class Blocked, class Nonblocked>
ATLAS_GLOBAL void kernel_copy_blocked_to_nonblocked_mdspan_perblock(const Blocked blocked, Nonblocked nonblocked) {
    idx_t npts   = nonblocked.extent(0);
    idx_t nproma = blocked.extent(blocked.rank()-1);
    idx_t nblk   = blocked.extent(0);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        idx_t work_units = nblk * nlev * nvar;
        auto traversal = KernelTraversal::current();
        for (idx_t unit = traversal.offset; unit < work_units; unit += traversal.stride) {
            idx_t jblk  = unit % nblk;
            idx_t entry = unit / nblk;
            idx_t jlev  = entry % nlev;
            idx_t jvar  = entry / nlev;
            idx_t base  = jblk * nproma;
            idx_t count = min(nproma, npts - base);
            for (idx_t jrof = 0; jrof < count; ++jrof) {
                nonblocked(base + jrof, jlev, jvar) = blocked(jblk, jvar, jlev, jrof);
            }
        }
    }
    else if constexpr (blocked.rank()==3) {
        idx_t nlev = nonblocked.extent(1);
        idx_t work_units = nblk * nlev;
        auto traversal = KernelTraversal::current();
        for (idx_t unit = traversal.offset; unit < work_units; unit += traversal.stride) {
            idx_t jblk = unit % nblk;
            idx_t jlev = unit / nblk;
            idx_t base = jblk * nproma;
            idx_t count = min(nproma, npts - base);
            for (idx_t jrof = 0; jrof < count; ++jrof) {
                nonblocked(base + jrof, jlev) = blocked(jblk, jlev, jrof);
            }
        }
    }
    else if constexpr (blocked.rank()==2) {
        auto traversal = KernelTraversal::current();
        for (idx_t jblk = traversal.offset; jblk < nblk; jblk += traversal.stride) {
            idx_t base = jblk * nproma;
            idx_t count = min(nproma, npts - base);
            for (idx_t jrof = 0; jrof < count; ++jrof) {
                nonblocked(base + jrof) = blocked(jblk, jrof);
            }
        }
    }
}

/// @brief Nonblocked-coalesced variant of the blocked-to-nonblocked kernel.
///
/// Maps consecutive threads to the innermost (unit-stride) dimension of the nonblocked target
/// instead of to the horizontal point index. This coalesces the nonblocked writes at the cost of
/// uncoalescing the blocked reads. For blocked->nonblocked the write side is the expensive side, so
/// trading coalesced reads for coalesced writes can be a net win. For rank 2 there is no interior
/// dimension, so this reduces to the flat mapping.
template <class Blocked, class Nonblocked>
ATLAS_GLOBAL void kernel_copy_blocked_to_nonblocked_mdspan_coalesced(const Blocked blocked, Nonblocked nonblocked) {
    idx_t npts   = nonblocked.extent(0);
    idx_t nproma = blocked.extent(blocked.rank()-1);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        idx_t work_size = npts * nlev * nvar;
        auto traversal = KernelTraversal::current();
        for (idx_t offset = traversal.offset; offset < work_size; offset += traversal.stride) {
            idx_t jvar  = offset % nvar;
            idx_t rem   = offset / nvar;
            idx_t jlev  = rem % nlev;
            idx_t point = rem / nlev;
            idx_t jblk  = point / nproma;
            idx_t jrof  = point - jblk * nproma;
            nonblocked(point, jlev, jvar) = blocked(jblk, jvar, jlev, jrof);
        }
    }
    else if constexpr (blocked.rank()==3) {
        idx_t nlev = nonblocked.extent(1);
        idx_t work_size = npts * nlev;
        auto traversal = KernelTraversal::current();
        for (idx_t offset = traversal.offset; offset < work_size; offset += traversal.stride) {
            idx_t jlev  = offset % nlev;
            idx_t point = offset / nlev;
            idx_t jblk  = point / nproma;
            idx_t jrof  = point - jblk * nproma;
            nonblocked(point, jlev) = blocked(jblk, jlev, jrof);
        }
    }
    else if constexpr (blocked.rank()==2) {
        auto traversal = KernelTraversal::current();
        for (idx_t point = traversal.offset; point < npts; point += traversal.stride) {
            idx_t jblk = point / nproma;
            idx_t jrof = point - jblk * nproma;
            nonblocked(point) = blocked(jblk, jrof);
        }
    }
}

template <class Nonblocked, class Blocked>
ATLAS_GLOBAL void kernel_copy_nonblocked_to_blocked_mdspan(Nonblocked nonblocked, Blocked blocked) {
    auto npts     = nonblocked.extent(0);
    auto nproma = blocked.extent(blocked.rank()-1);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        idx_t work_size = npts * nlev * nvar;
        auto traversal = KernelTraversal::current();
        for (idx_t offset = traversal.offset; offset < work_size; offset += traversal.stride) {
            idx_t point = offset % npts;
            idx_t entry = offset / npts;
            idx_t jlev  = entry % nlev;
            idx_t jvar  = entry / nlev;
            idx_t jblk  = point / nproma;
            idx_t jrof  = point - jblk * nproma;
            blocked(jblk, jvar, jlev, jrof) = nonblocked(point, jlev, jvar);
        }
    }
    else if constexpr (blocked.rank()==3) {
        idx_t nlev = nonblocked.extent(1);
        idx_t work_size = npts * nlev;
        auto traversal = KernelTraversal::current();
        for (idx_t offset = traversal.offset; offset < work_size; offset += traversal.stride) {
            idx_t point = offset % npts;
            idx_t jlev  = offset / npts;
            idx_t jblk  = point / nproma;
            idx_t jrof  = point - jblk * nproma;
            blocked(jblk, jlev, jrof) = nonblocked(point, jlev);
        }
    }
    else if constexpr (blocked.rank()==2) {
        auto traversal = KernelTraversal::current();
        for (idx_t point = traversal.offset; point < npts; point += traversal.stride) {
            idx_t jblk = point / nproma;
            idx_t jrof = point - jblk * nproma;
            blocked(jblk, jrof) = nonblocked(point);
        }
    }
}

/// @brief Nonblocked-coalesced variant of the nonblocked-to-blocked kernel.
///
/// Mirror of kernel_copy_blocked_to_nonblocked_mdspan_coalesced: maps consecutive threads to the
/// innermost (unit-stride) dimension of the nonblocked source instead of to the horizontal point
/// index. This coalesces the nonblocked reads at the cost of uncoalescing the blocked writes, but
/// the resulting blocked writes stride by nproma (a small stride) rather than by nlev, so the
/// scattered side is much tighter than in the flat nproma_outermost kernel. For rank 2 there is no
/// interior dimension, so this reduces to the flat mapping.
template <class Nonblocked, class Blocked>
ATLAS_GLOBAL void kernel_copy_nonblocked_to_blocked_mdspan_coalesced(Nonblocked nonblocked, Blocked blocked) {
    idx_t npts   = nonblocked.extent(0);
    idx_t nproma = blocked.extent(blocked.rank()-1);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        idx_t work_size = npts * nlev * nvar;
        auto traversal = KernelTraversal::current();
        for (idx_t offset = traversal.offset; offset < work_size; offset += traversal.stride) {
            idx_t jvar  = offset % nvar;
            idx_t rem   = offset / nvar;
            idx_t jlev  = rem % nlev;
            idx_t point = rem / nlev;
            idx_t jblk  = point / nproma;
            idx_t jrof  = point - jblk * nproma;
            blocked(jblk, jvar, jlev, jrof) = nonblocked(point, jlev, jvar);
        }
    }
    else if constexpr (blocked.rank()==3) {
        idx_t nlev = nonblocked.extent(1);
        idx_t work_size = npts * nlev;
        auto traversal = KernelTraversal::current();
        for (idx_t offset = traversal.offset; offset < work_size; offset += traversal.stride) {
            idx_t jlev  = offset % nlev;
            idx_t point = offset / nlev;
            idx_t jblk  = point / nproma;
            idx_t jrof  = point - jblk * nproma;
            blocked(jblk, jlev, jrof) = nonblocked(point, jlev);
        }
    }
    else if constexpr (blocked.rank()==2) {
        auto traversal = KernelTraversal::current();
        for (idx_t point = traversal.offset; point < npts; point += traversal.stride) {
            idx_t jblk = point / nproma;
            idx_t jrof = point - jblk * nproma;
            blocked(jblk, jrof) = nonblocked(point);
        }
    }
}

/// @brief Per-block (nproma_innermost) variant of the nonblocked-to-blocked kernel.
///
/// Inverse of kernel_copy_blocked_to_nonblocked_mdspan_perblock: each work item owns one blocked
/// chunk (jblk) for a fixed (jlev, jvar) and loops over the contiguous nproma dimension innermost,
/// giving unit-stride writes on the blocked side and hoisting the per-element div/mod.
template <class Nonblocked, class Blocked>
ATLAS_GLOBAL void kernel_copy_nonblocked_to_blocked_mdspan_perblock(Nonblocked nonblocked, Blocked blocked) {
    idx_t npts   = nonblocked.extent(0);
    idx_t nproma = blocked.extent(blocked.rank()-1);
    idx_t nblk   = blocked.extent(0);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        idx_t work_units = nblk * nlev * nvar;
        auto traversal = KernelTraversal::current();
        for (idx_t unit = traversal.offset; unit < work_units; unit += traversal.stride) {
            idx_t jblk  = unit % nblk;
            idx_t entry = unit / nblk;
            idx_t jlev  = entry % nlev;
            idx_t jvar  = entry / nlev;
            idx_t base  = jblk * nproma;
            idx_t count = min(nproma, npts - base);
            for (idx_t jrof = 0; jrof < count; ++jrof) {
                blocked(jblk, jvar, jlev, jrof) = nonblocked(base + jrof, jlev, jvar);
            }
        }
    }
    else if constexpr (blocked.rank()==3) {
        idx_t nlev = nonblocked.extent(1);
        idx_t work_units = nblk * nlev;
        auto traversal = KernelTraversal::current();
        for (idx_t unit = traversal.offset; unit < work_units; unit += traversal.stride) {
            idx_t jblk = unit % nblk;
            idx_t jlev = unit / nblk;
            idx_t base = jblk * nproma;
            idx_t count = min(nproma, npts - base);
            for (idx_t jrof = 0; jrof < count; ++jrof) {
                blocked(jblk, jlev, jrof) = nonblocked(base + jrof, jlev);
            }
        }
    }
    else if constexpr (blocked.rank()==2) {
        auto traversal = KernelTraversal::current();
        for (idx_t jblk = traversal.offset; jblk < nblk; jblk += traversal.stride) {
            idx_t base = jblk * nproma;
            idx_t count = min(nproma, npts - base);
            for (idx_t jrof = 0; jrof < count; ++jrof) {
                blocked(jblk, jrof) = nonblocked(base + jrof);
            }
        }
    }
}

template <class Nonblocked, class Blocked>
/**
 * @brief Launch or run the device copy from a nonblocked view to a blocked view.
 *
 * @param nonblocked Source device-accessible Atlas view (`atlas::View`/`atlas::ArrayView`) or
 *        mdspan-like view with rank one less than `blocked`.
 * @param blocked Target device-accessible Atlas view (`atlas::View`/`atlas::ArrayView`) or
 *        mdspan-like view with rank 2, 3, or 4.
 *
 * @pre `nonblocked.rank() == blocked.rank() - 1`.
 * @pre Shared dimensions must match; for rank 4, `nonblocked.extent(1) == blocked.extent(2)`
 *      and `nonblocked.extent(2) == blocked.extent(1)`.
 */
void device_copy_nonblocked_to_blocked_impl(Nonblocked nonblocked, Blocked blocked) {
    if constexpr(blocked.rank()==4) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
    }
    else if constexpr (blocked.rank()==3) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
    }
    const RelayoutLoopOrder order = resolve_loop_order(relayout_loop_order(), CopyDirection::nonblocked_to_blocked);
    #if HIC_COMPILER
    if (order == RelayoutLoopOrder::nproma_innermost) {
        idx_t work_units = blocked.extent(0);
        if constexpr(blocked.rank()==4) {
            work_units *= nonblocked.extent(1) * nonblocked.extent(2);
        }
        else if constexpr(blocked.rank()==3) {
            work_units *= nonblocked.extent(1);
        }
        kernel_copy_nonblocked_to_blocked_mdspan_perblock<<<blocks_per_grid(work_units),threads_per_block()>>>(nonblocked, blocked);
    }
    else if (order == RelayoutLoopOrder::nonblocked_coalesced) {
        idx_t work_size = nonblocked.extent(0);
        if constexpr(blocked.rank()==4) {
            work_size *= nonblocked.extent(1) * nonblocked.extent(2);
        }
        else if constexpr(blocked.rank()==3) {
            work_size *= nonblocked.extent(1);
        }
        kernel_copy_nonblocked_to_blocked_mdspan_coalesced<<<blocks_per_grid(work_size),threads_per_block()>>>(nonblocked, blocked);
    }
    else {
        idx_t work_size = nonblocked.extent(0);
        if constexpr(blocked.rank()==4) {
            work_size *= nonblocked.extent(1) * nonblocked.extent(2);
        }
        else if constexpr(blocked.rank()==3) {
            work_size *= nonblocked.extent(1);
        }
        kernel_copy_nonblocked_to_blocked_mdspan<<<blocks_per_grid(work_size),threads_per_block()>>>(nonblocked, blocked);
    }
    HIC_CHECK_KERNEL_LAUNCH();
    #else
    if (order == RelayoutLoopOrder::nproma_innermost) {
        kernel_copy_nonblocked_to_blocked_mdspan_perblock(nonblocked, blocked);
    }
    else if (order == RelayoutLoopOrder::nonblocked_coalesced) {
        kernel_copy_nonblocked_to_blocked_mdspan_coalesced(nonblocked, blocked);
    }
    else {
        kernel_copy_nonblocked_to_blocked_mdspan(nonblocked, blocked);
    }
    #endif
}

template <class Blocked, class Nonblocked>
/**
 * @brief Launch or run the device copy from a blocked view to a nonblocked view.
 *
 * @param blocked Source device-accessible Atlas view (`atlas::View`/`atlas::ArrayView`) or
 *        mdspan-like view with rank 2, 3, or 4.
 * @param nonblocked Target device-accessible Atlas view (`atlas::View`/`atlas::ArrayView`) or
 *        mdspan-like view with rank one less than `blocked`.
 *
 * @pre `nonblocked.rank() == blocked.rank() - 1`.
 * @pre Shared dimensions must match; for rank 4, `nonblocked.extent(1) == blocked.extent(2)`
 *      and `nonblocked.extent(2) == blocked.extent(1)`.
 */
void device_copy_blocked_to_nonblocked_impl(Blocked blocked, Nonblocked nonblocked) {
    ATLAS_ASSERT(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
    }
    else if constexpr (blocked.rank()==3) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
    }
    else if constexpr (blocked.rank()==2) {
    }
    else {
        ATLAS_THROW_EXCEPTION("transposition not implemented");
    }
    const RelayoutLoopOrder order = resolve_loop_order(relayout_loop_order(), CopyDirection::blocked_to_nonblocked);
    #if HIC_COMPILER
    if (order == RelayoutLoopOrder::nproma_innermost) {
        idx_t work_units = blocked.extent(0);
        if constexpr(blocked.rank()==4) {
            work_units *= nonblocked.extent(1) * nonblocked.extent(2);
        }
        else if constexpr(blocked.rank()==3) {
            work_units *= nonblocked.extent(1);
        }
        kernel_copy_blocked_to_nonblocked_mdspan_perblock<<<blocks_per_grid(work_units),threads_per_block()>>>(blocked, nonblocked);
    }
    else if (order == RelayoutLoopOrder::nonblocked_coalesced) {
        idx_t work_size = nonblocked.extent(0);
        if constexpr(blocked.rank()==4) {
            work_size *= nonblocked.extent(1) * nonblocked.extent(2);
        }
        else if constexpr(blocked.rank()==3) {
            work_size *= nonblocked.extent(1);
        }
        kernel_copy_blocked_to_nonblocked_mdspan_coalesced<<<blocks_per_grid(work_size),threads_per_block()>>>(blocked, nonblocked);
    }
    else {
        idx_t work_size = nonblocked.extent(0);
        if constexpr(blocked.rank()==4) {
            work_size *= nonblocked.extent(1) * nonblocked.extent(2);
        }
        else if constexpr(blocked.rank()==3) {
            work_size *= nonblocked.extent(1);
        }
        kernel_copy_blocked_to_nonblocked_mdspan<<<blocks_per_grid(work_size),threads_per_block()>>>(blocked, nonblocked);
    }
    HIC_CHECK_KERNEL_LAUNCH();
    #else
    if (order == RelayoutLoopOrder::nproma_innermost) {
        kernel_copy_blocked_to_nonblocked_mdspan_perblock(blocked, nonblocked);
    }
    else if (order == RelayoutLoopOrder::nonblocked_coalesced) {
        kernel_copy_blocked_to_nonblocked_mdspan_coalesced(blocked, nonblocked);
    }
    else {
        kernel_copy_blocked_to_nonblocked_mdspan(blocked, nonblocked);
    }
    #endif
}

template <class BlockedIn, class BlockedOut>
ATLAS_GLOBAL void kernel_copy_blocked_to_blocked_mdspan(BlockedIn blocked_in, BlockedOut blocked_out) {
    auto nblks_in  = blocked_in.extent(0);
    auto nproma_in = blocked_in.extent(blocked_in.rank()-1);
    auto nblks_out  = blocked_out.extent(0);
    auto nproma_out = blocked_out.extent(blocked_out.rank()-1);
    auto total_points_in  = nblks_in * nproma_in;
    auto total_points_out = nblks_out * nproma_out;
    auto total_points     = min(total_points_in, total_points_out);

    if constexpr (blocked_in.rank()==4) {
        idx_t nlev = blocked_in.extent(1);
        idx_t nvar = blocked_in.extent(2);
        idx_t work_size = total_points * nlev * nvar;
        auto traversal = KernelTraversal::current();
        for (idx_t offset = traversal.offset; offset < work_size; offset += traversal.stride) {
            idx_t point = offset % total_points;
            idx_t entry = offset / total_points;
            idx_t jlev  = entry % nlev;
            idx_t jvar  = entry / nlev;
            idx_t jblk_in   = point / nproma_in;
            idx_t jrof_in   = point - jblk_in * nproma_in;
            idx_t jblk_out  = point / nproma_out;
            idx_t jrof_out  = point - jblk_out * nproma_out;
            blocked_out(jblk_out, jlev, jvar, jrof_out) = blocked_in(jblk_in, jlev, jvar, jrof_in);
        }
    }
    else if constexpr (blocked_in.rank()==3) {
        idx_t nlev = blocked_in.extent(1);
        idx_t work_size = total_points * nlev;
        auto traversal = KernelTraversal::current();
        for (idx_t offset = traversal.offset; offset < work_size; offset += traversal.stride) {
            idx_t point = offset % total_points;
            idx_t jlev  = offset / total_points;
            idx_t jblk_in   = point / nproma_in;
            idx_t jrof_in   = point - jblk_in * nproma_in;
            idx_t jblk_out  = point / nproma_out;
            idx_t jrof_out  = point - jblk_out * nproma_out;
            blocked_out(jblk_out, jlev, jrof_out) = blocked_in(jblk_in, jlev, jrof_in);
        }
    }
    else if constexpr (blocked_in.rank()==2) {
        auto traversal = KernelTraversal::current();
        for (idx_t point = traversal.offset; point < total_points; point += traversal.stride) {
            idx_t jblk_in   = point / nproma_in;
            idx_t jrof_in   = point - jblk_in * nproma_in;
            idx_t jblk_out  = point / nproma_out;
            idx_t jrof_out  = point - jblk_out * nproma_out;
            blocked_out(jblk_out, jrof_out) = blocked_in(jblk_in, jrof_in);
        }
    }
}

/// @brief Per-block (nproma_innermost) variant of the blocked-to-blocked kernel.
///
/// Each work item owns one output blocked chunk (jblk_out) for a fixed (jlev, jvar) and loops over
/// the contiguous output nproma dimension innermost, giving unit-stride writes on the output side.
/// The input indices are still computed per element because the two fields may use a different
/// nproma, so the input side remains a general gather.
template <class BlockedIn, class BlockedOut>
ATLAS_GLOBAL void kernel_copy_blocked_to_blocked_mdspan_perblock(BlockedIn blocked_in, BlockedOut blocked_out) {
    idx_t nproma_in  = blocked_in.extent(blocked_in.rank()-1);
    idx_t nproma_out = blocked_out.extent(blocked_out.rank()-1);
    idx_t nblks_out  = blocked_out.extent(0);
    idx_t total_points_in  = blocked_in.extent(0) * nproma_in;
    idx_t total_points_out = nblks_out * nproma_out;
    idx_t total_points     = min(total_points_in, total_points_out);
    if constexpr (blocked_in.rank()==4) {
        idx_t nlev = blocked_in.extent(1);
        idx_t nvar = blocked_in.extent(2);
        idx_t work_units = nblks_out * nlev * nvar;
        auto traversal = KernelTraversal::current();
        for (idx_t unit = traversal.offset; unit < work_units; unit += traversal.stride) {
            idx_t jblk_out = unit % nblks_out;
            idx_t entry    = unit / nblks_out;
            idx_t jlev     = entry % nlev;
            idx_t jvar     = entry / nlev;
            idx_t base_out = jblk_out * nproma_out;
            for (idx_t jrof_out = 0; jrof_out < nproma_out; ++jrof_out) {
                idx_t point = base_out + jrof_out;
                if (point >= total_points) break;
                idx_t jblk_in = point / nproma_in;
                idx_t jrof_in = point - jblk_in * nproma_in;
                blocked_out(jblk_out, jlev, jvar, jrof_out) = blocked_in(jblk_in, jlev, jvar, jrof_in);
            }
        }
    }
    else if constexpr (blocked_in.rank()==3) {
        idx_t nlev = blocked_in.extent(1);
        idx_t work_units = nblks_out * nlev;
        auto traversal = KernelTraversal::current();
        for (idx_t unit = traversal.offset; unit < work_units; unit += traversal.stride) {
            idx_t jblk_out = unit % nblks_out;
            idx_t jlev     = unit / nblks_out;
            idx_t base_out = jblk_out * nproma_out;
            for (idx_t jrof_out = 0; jrof_out < nproma_out; ++jrof_out) {
                idx_t point = base_out + jrof_out;
                if (point >= total_points) break;
                idx_t jblk_in = point / nproma_in;
                idx_t jrof_in = point - jblk_in * nproma_in;
                blocked_out(jblk_out, jlev, jrof_out) = blocked_in(jblk_in, jlev, jrof_in);
            }
        }
    }
    else if constexpr (blocked_in.rank()==2) {
        auto traversal = KernelTraversal::current();
        for (idx_t jblk_out = traversal.offset; jblk_out < nblks_out; jblk_out += traversal.stride) {
            idx_t base_out = jblk_out * nproma_out;
            for (idx_t jrof_out = 0; jrof_out < nproma_out; ++jrof_out) {
                idx_t point = base_out + jrof_out;
                if (point >= total_points) break;
                idx_t jblk_in = point / nproma_in;
                idx_t jrof_in = point - jblk_in * nproma_in;
                blocked_out(jblk_out, jrof_out) = blocked_in(jblk_in, jrof_in);
            }
        }
    }
}

template <class BlockedIn, class BlockedOut>
/**
 * @brief Launch or run the device copy between two blocked views.
 *
 * @param blocked_in Source device-accessible Atlas view (`atlas::View`/`atlas::ArrayView`) or
 *        mdspan-like blocked view with rank 2, 3, or 4.
 * @param blocked_out Target device-accessible Atlas view (`atlas::View`/`atlas::ArrayView`) or
 *        mdspan-like blocked view with the same rank and value type as `blocked_in`.
 *
 * @pre `blocked_in.rank() == blocked_out.rank()`.
 * @pre Non-horizontal dimensions must match.
 */
void device_copy_blocked_to_blocked_impl(BlockedIn blocked_in, BlockedOut blocked_out) {
    static_assert(blocked_in.rank() == blocked_out.rank());
    if constexpr (blocked_in.rank()==4) {
        ATLAS_ASSERT(blocked_in.extent(1) == blocked_out.extent(1));
        ATLAS_ASSERT(blocked_in.extent(2) == blocked_out.extent(2));
    }
    else
    if constexpr (blocked_in.rank()==3) {
        ATLAS_ASSERT(blocked_in.extent(1) == blocked_out.extent(1));
    }
    #if HIC_COMPILER
    if (relayout_loop_order() == RelayoutLoopOrder::nproma_innermost) {
        idx_t work_units = blocked_out.extent(0);
        if constexpr(blocked_in.rank()==4) {
            work_units *= blocked_in.extent(1) * blocked_in.extent(2);
        }
        else if constexpr(blocked_in.rank()==3) {
            work_units *= blocked_in.extent(1);
        }
        kernel_copy_blocked_to_blocked_mdspan_perblock<<<blocks_per_grid(work_units),threads_per_block()>>>(blocked_in, blocked_out);
    }
    else {
        idx_t input_points = blocked_in.extent(0) * blocked_in.extent(blocked_in.rank()-1);
        idx_t output_points = blocked_out.extent(0) * blocked_out.extent(blocked_out.rank()-1);
        idx_t work_size = input_points < output_points ? input_points : output_points;
        if constexpr(blocked_in.rank()==4) {
            work_size *= blocked_in.extent(1) * blocked_in.extent(2);
        }
        else if constexpr(blocked_in.rank()==3) {
            work_size *= blocked_in.extent(1);
        }
        kernel_copy_blocked_to_blocked_mdspan<<<blocks_per_grid(work_size),threads_per_block()>>>(blocked_in, blocked_out);
    }
    HIC_CHECK_KERNEL_LAUNCH();
    #else
    if (relayout_loop_order() == RelayoutLoopOrder::nproma_innermost) {
        kernel_copy_blocked_to_blocked_mdspan_perblock(blocked_in, blocked_out);
    }
    else {
        kernel_copy_blocked_to_blocked_mdspan(blocked_in, blocked_out);
    }
    #endif
}

template <class Nonblocked, class Blocked>
void device_copy_nonblocked_to_blocked_mdspan(Nonblocked nonblocked, Blocked blocked) {
    if (relayout_implementation() == RelayoutImplementation::mdspan) {
        return dispatch_relayout_mdspan<ViewType::Nonblocked, ViewType::Blocked>(nonblocked, blocked, [&](const auto& nonblocked_view, const auto& blocked_view) {
            return device_copy_nonblocked_to_blocked_impl(nonblocked_view, blocked_view);
        });
    }
    return device_copy_nonblocked_to_blocked_impl(nonblocked, blocked);
}

template <class Blocked, class Nonblocked>
void device_copy_blocked_to_nonblocked_mdspan(Blocked blocked, Nonblocked nonblocked) {
    if (relayout_implementation() == RelayoutImplementation::mdspan) {
        return dispatch_relayout_mdspan<ViewType::Blocked, ViewType::Nonblocked>(blocked, nonblocked, [&](const auto& blocked_view, const auto& nonblocked_view) {
            return device_copy_blocked_to_nonblocked_impl(blocked_view, nonblocked_view);
        });
    }
    return device_copy_blocked_to_nonblocked_impl(blocked, nonblocked);
}

template <class BlockedIn, class BlockedOut>
void device_copy_blocked_to_blocked_mdspan(BlockedIn blocked_in, BlockedOut blocked_out) {
    if (relayout_implementation() == RelayoutImplementation::mdspan) {
        return dispatch_relayout_mdspan<ViewType::Blocked, ViewType::Blocked>(blocked_in, blocked_out, [&](const auto& blocked_in_view, const auto& blocked_out_view) {
            return device_copy_blocked_to_blocked_impl(blocked_in_view, blocked_out_view);
        });
    }
    return device_copy_blocked_to_blocked_impl(blocked_in, blocked_out);
}

} //namespace atlas


//-----------------------------------------------------------------------
// Explicit instantiation
//
// The actual instantiations live in the per-type translation units
// relayout_on_device_{int,float,long,double}.hic, which include this file and
// expand ATLAS_RELAYOUT_DEVICE_EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK for a
// single type. This mirrors the host split (relayout_on_host_{int,...}.cc).

#define ATLAS_RELAYOUT_DEVICE_EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(TYPE, BLOCKED_RANK) \
    template void atlas::device_copy_blocked_to_nonblocked_mdspan<atlas::array::ArrayView<const TYPE,BLOCKED_RANK>,atlas::array::ArrayView<TYPE,BLOCKED_RANK-1>>(atlas::array::ArrayView<const TYPE,BLOCKED_RANK>, atlas::array::ArrayView<TYPE,BLOCKED_RANK-1>); \
    template void atlas::device_copy_nonblocked_to_blocked_mdspan<atlas::array::ArrayView<const TYPE,BLOCKED_RANK-1>,atlas::array::ArrayView<TYPE,BLOCKED_RANK>>(atlas::array::ArrayView<const TYPE,BLOCKED_RANK-1>, atlas::array::ArrayView<TYPE,BLOCKED_RANK>); \
    template void atlas::device_copy_blocked_to_blocked_mdspan<atlas::array::ArrayView<const TYPE,BLOCKED_RANK>,atlas::array::ArrayView<TYPE,BLOCKED_RANK>>(atlas::array::ArrayView<const TYPE,BLOCKED_RANK>, atlas::array::ArrayView<TYPE,BLOCKED_RANK>);
