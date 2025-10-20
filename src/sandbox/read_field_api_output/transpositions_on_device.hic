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
#include "atlas/array.h"
#include "hic/hic.h"

#include "atlas/runtime/Log.h"

namespace atlas {

#define USE_MDSPAN 0

ATLAS_DEVICE inline idx_t min(idx_t a, idx_t b) {
    return a < b ? a : b;
}

template <class Blocked, class Nonblocked>
ATLAS_GLOBAL void kernel_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    auto np     = nonblocked.extent(0);
    auto nblks  = blocked.extent(0);
    auto nproma = blocked.extent(blocked.rank()-1);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            auto nrof = min(np - jpbegin, nproma);
            for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        auto jp = jpbegin+jrof;
                        nonblocked(jp, jlev, jvar) = blocked(jblk, jvar, jlev, jrof);
                    }
                }
            }
        }
    }
    else if constexpr (blocked.rank()==3) {
        idx_t nlev = nonblocked.extent(1);
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            auto nrof = min(np - jpbegin, nproma);
            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    auto jp = jpbegin+jrof;
                    nonblocked(jp, jlev) = blocked(jblk, jlev, jrof);
                }
            }
        }
    }
    else if constexpr (blocked.rank()==2) {
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            auto nrof = min(np - jpbegin, nproma);
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                auto jp = jpbegin+jrof;
                nonblocked(jp) = blocked(jblk, jrof);
            }
        }
    }
}

template <class Nonblocked, class Blocked>
ATLAS_GLOBAL void kernel_copy_nonblocked_to_blocked_mdspan(Nonblocked nonblocked, Blocked blocked) {
    auto np     = nonblocked.extent(0);
    auto nblks  = blocked.extent(0);
    auto nproma = blocked.extent(blocked.rank()-1);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            auto nrof = min(np - jpbegin, nproma);
            for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        auto jp = jpbegin+jrof;
                        blocked(jblk, jvar, jlev, jrof) = nonblocked(jp, jlev, jvar);
                    }
                }
            }
        }
    }
    else if constexpr (blocked.rank()==3) {
        idx_t nlev = nonblocked.extent(1);
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            auto nrof = min(np - jpbegin, nproma);
            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                    auto jp = jpbegin+jrof;
                     blocked(jblk, jlev, jrof) = nonblocked(jp, jlev);
                }
            }
        }
    }
    else if constexpr (blocked.rank()==2) {
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            auto nrof = min(np - jpbegin, nproma);
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                auto jp = jpbegin+jrof;
                blocked(jblk, jrof) = nonblocked(jp);
            }
        }
    }
}

template <class Nonblocked, class Blocked>
void device_copy_nonblocked_to_blocked_mdspan(Nonblocked nonblocked, Blocked blocked) {
    if constexpr(blocked.rank()==4) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
    }
    else if constexpr (blocked.rank()==3) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
    }
    #if HIC_COMPILER
    HIC_CALL( kernel_copy_nonblocked_to_blocked_mdspan<<<1,1>>>(nonblocked, blocked) );
    #else
    kernel_copy_nonblocked_to_blocked_mdspan(nonblocked, blocked);
    #endif
}

template <class Blocked, class Nonblocked>
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
    HIC_CALL( kernel_copy_blocked_to_nonblocked_mdspan<<<1,1>>>(blocked, nonblocked) );
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
    if constexpr (blocked_in.rank()==3) {
        idx_t nlev = blocked_in.extent(1);
        idx_t jp_out = 0;
        for (idx_t jblk_in = 0; jblk_in < nblks_in; ++jblk_in) {
            for (idx_t jrof_in = 0; jrof_in < nproma_in; ++jrof_in, ++jp_out) {
                idx_t jblk_out = jp_out / nproma_out;
                if (jblk_out < nblks_out) {
                    idx_t jrof_out = jp_out - jblk_out * nproma_out;
                    if (jrof_out < nproma_out) {
                        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                            blocked_out(jblk_out, jlev, jrof_out) = blocked_in(jblk_in, jlev, jrof_in);
                        }
                    }
                }
            }
        }
    }
    else if constexpr (blocked_in.rank()==2) {
        idx_t jp_out = 0;
        for (idx_t jblk_in = 0; jblk_in < nblks_in; ++jblk_in) {
            for (idx_t jrof_in = 0; jrof_in < nproma_in; ++jrof_in, ++jp_out) {
                idx_t jblk_out = jp_out / nproma_out;
                if (jblk_out < nblks_out) {
                    idx_t jrof_out = jp_out - jblk_out * nproma_out;
                    if (jrof_out < nproma_out) {
                        blocked_out(jblk_out, jrof_out) = blocked_in(jblk_in, jrof_in);
                    }
                }
            }
        }
    }
}

template <class BlockedIn, class BlockedOut>
void device_copy_blocked_to_blocked_mdspan(BlockedIn blocked_in, BlockedOut blocked_out) {
    static_assert(blocked_in.rank() == blocked_out.rank());
    if constexpr (blocked_in.rank()==3) {
        ATLAS_ASSERT(blocked_in.extent(1) == blocked_out.extent(1));
    }
    #if HIC_COMPILER
    HIC_CALL(kernel_copy_blocked_to_blocked_mdspan<<<1,1>>>(blocked_in, blocked_out));
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
    template void device_copy_nonblocked_to_blocked_mdspan<array::ArrayView<const TYPE,BLOCKED_RANK-1>,array::ArrayView<TYPE,BLOCKED_RANK>>(array::ArrayView<const TYPE,BLOCKED_RANK-1>, array::ArrayView<TYPE,BLOCKED_RANK>); \
    template void device_copy_blocked_to_blocked_mdspan<array::ArrayView<const TYPE,BLOCKED_RANK>,array::ArrayView<TYPE,BLOCKED_RANK>>(array::ArrayView<const TYPE,BLOCKED_RANK>, array::ArrayView<TYPE,BLOCKED_RANK>);

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
