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

/**
 * @file relayout_on_device.hic
 * @brief Device implementations for copying Atlas data between blocked and nonblocked layouts.
 *
 * Device relayout expects device-accessible Atlas views (`atlas::View`/`atlas::ArrayView`) or
 * mdspan-like views with the rank, datatype, and shape requirements documented in
 * `atlas/util/relayout.h`.
 */

#include "hic/hic.h"
#include "atlas/array.h"
#include "atlas/util/relayout.h"

#include "atlas/runtime/Log.h"

namespace atlas {

#define USE_MDSPAN 0

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
void device_copy_nonblocked_to_blocked_mdspan(Nonblocked nonblocked, Blocked blocked) {
    if constexpr(blocked.rank()==4) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
    }
    else if constexpr (blocked.rank()==3) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
    }
    #if HIC_COMPILER
    idx_t work_size = nonblocked.extent(0);
    if constexpr(blocked.rank()==4) {
        work_size *= nonblocked.extent(1) * nonblocked.extent(2);
    }
    else if constexpr(blocked.rank()==3) {
        work_size *= nonblocked.extent(1);
    }
    kernel_copy_nonblocked_to_blocked_mdspan<<<blocks_per_grid(work_size),threads_per_block()>>>(nonblocked, blocked);
    HIC_CHECK_KERNEL_LAUNCH();
    #else
    kernel_copy_nonblocked_to_blocked_mdspan(nonblocked, blocked);
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
void device_copy_blocked_to_nonblocked_mdspan(Blocked blocked, Nonblocked nonblocked) {
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
    #if HIC_COMPILER
    idx_t work_size = nonblocked.extent(0);
    if constexpr(blocked.rank()==4) {
        work_size *= nonblocked.extent(1) * nonblocked.extent(2);
    }
    else if constexpr(blocked.rank()==3) {
        work_size *= nonblocked.extent(1);
    }
    kernel_copy_blocked_to_nonblocked_mdspan<<<blocks_per_grid(work_size),threads_per_block()>>>(blocked, nonblocked);
    HIC_CHECK_KERNEL_LAUNCH();
    #else
    kernel_copy_blocked_to_nonblocked_mdspan(blocked, nonblocked);
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
void device_copy_blocked_to_blocked_mdspan(BlockedIn blocked_in, BlockedOut blocked_out) {
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
    HIC_CHECK_KERNEL_LAUNCH();
    #else
    kernel_copy_blocked_to_blocked_mdspan(blocked_in, blocked_out);
    #endif

}

} //namespace atlas


//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {

#define EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(TYPE, BLOCKED_RANK) \
    template void device_copy_blocked_to_nonblocked_mdspan<array::ArrayView<const TYPE,BLOCKED_RANK>,array::ArrayView<TYPE,BLOCKED_RANK-1>>(array::ArrayView<const TYPE,BLOCKED_RANK>, array::ArrayView<TYPE,BLOCKED_RANK-1>); \
    template void device_copy_blocked_to_nonblocked_mdspan<mdspan<const TYPE,dims<BLOCKED_RANK>,layout_stride>,mdspan<TYPE,dims<BLOCKED_RANK-1>,layout_stride>>(mdspan<const TYPE,dims<BLOCKED_RANK>,layout_stride>, mdspan<TYPE,dims<BLOCKED_RANK-1>,layout_stride>); \
\
    template void device_copy_nonblocked_to_blocked_mdspan<array::ArrayView<const TYPE,BLOCKED_RANK-1>,array::ArrayView<TYPE,BLOCKED_RANK>>(array::ArrayView<const TYPE,BLOCKED_RANK-1>, array::ArrayView<TYPE,BLOCKED_RANK>); \
    template void device_copy_nonblocked_to_blocked_mdspan<mdspan<const TYPE,dims<BLOCKED_RANK-1>,layout_stride>,mdspan<TYPE,dims<BLOCKED_RANK>,layout_stride>>(mdspan<const TYPE,dims<BLOCKED_RANK-1>,layout_stride>, mdspan<TYPE,dims<BLOCKED_RANK>,layout_stride>); \
\
    template void device_copy_blocked_to_blocked_mdspan<array::ArrayView<const TYPE,BLOCKED_RANK>,array::ArrayView<TYPE,BLOCKED_RANK>>(array::ArrayView<const TYPE,BLOCKED_RANK>, array::ArrayView<TYPE,BLOCKED_RANK>); \
    template void device_copy_blocked_to_blocked_mdspan<mdspan<const TYPE,dims<BLOCKED_RANK>,layout_stride>,mdspan<TYPE,dims<BLOCKED_RANK>,layout_stride>>(mdspan<const TYPE,dims<BLOCKED_RANK>,layout_stride>, mdspan<TYPE,dims<BLOCKED_RANK>,layout_stride>);

#define EXPLICIT_TEMPLATE_INSTATIATION(RANK)                \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(double, RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(float , RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(int   , RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(long  , RANK)


EXPLICIT_TEMPLATE_INSTATIATION(2)
EXPLICIT_TEMPLATE_INSTATIATION(3)
EXPLICIT_TEMPLATE_INSTATIATION(4)

#undef EXPLICIT_TEMPLATE_INSTATIATION_TYPE_RANK
#undef EXPLICIT_TEMPLATE_INSTATIATION

}  // namespace atlas
