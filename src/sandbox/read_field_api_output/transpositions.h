/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once
#include "atlas/array.h"
#include "atlas/field.h"

#include "atlas/functionspace/BlockStructuredColumns.h"

namespace atlas {

#define USE_MDSPAN 0 // Using mdspan is slower on macos!

template <class Nonblocked, class Blocked>
void host_copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked) {
    auto np     = nonblocked.extent(0);
    auto nblks  = blocked.extent(0);
    auto nproma = blocked.extent(blocked.rank()-1);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            auto nrof = std::min(np - jpbegin, nproma);
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
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
        idx_t nlev = nonblocked.extent(1);
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            auto nrof = std::min(np - jpbegin, nproma);
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
            auto nrof = std::min(np - jpbegin, nproma);
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                auto jp = jpbegin+jrof;
                blocked(jblk, jrof) = nonblocked(jp);
            }
        }
    }
    else {
        ATLAS_THROW_EXCEPTION("transposition not implemented");
    }
}

template <class Blocked, class Nonblocked>
void host_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    auto np     = nonblocked.extent(0);
    auto nblks  = blocked.extent(0);
    auto nproma = blocked.extent(blocked.rank()-1);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            auto nrof = std::min(np - jpbegin, nproma);
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
        ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
        idx_t nlev = nonblocked.extent(1);
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            auto nrof = std::min(np - jpbegin, nproma);
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
            auto nrof = std::min(np - jpbegin, nproma);
            for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                auto jp = jpbegin+jrof;
                nonblocked(jp) = blocked(jblk, jrof);
            }
        }
    }
    else {
        ATLAS_THROW_EXCEPTION("transposition not implemented");
    }
}

template <class BlockedIn, class BlockedOut>
void host_copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out) {
    auto nblks_in  = blocked_in.extent(0);
    auto nproma_in = blocked_in.extent(blocked_in.rank()-1);
    auto nblks_out  = blocked_out.extent(0);
    auto nproma_out = blocked_out.extent(blocked_out.rank()-1);
    static_assert(blocked_in.rank() == blocked_out.rank());
    if constexpr (blocked_in.rank()==3) {
        ATLAS_ASSERT(blocked_in.extent(1) == blocked_out.extent(1));
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
    else {
        ATLAS_THROW_EXCEPTION("transposition not implemented for rank " << blocked_in.rank());
    }
}

#if ATLAS_HAVE_GPU
template <class Nonblocked, class Blocked>
void device_copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked) {
    if (pluto::devices() == 0) {
        return host_copy_nonblocked_to_blocked_mdspan(nonblocked, blocked);
    }
    ATLAS_NOTIMPLEMENTED;
}

template <class Blocked, class Nonblocked>
void device_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    if (pluto::devices() == 0) {
        return host_copy_blocked_to_nonblocked_mdspan(blocked, nonblocked);
    }
    ATLAS_NOTIMPLEMENTED;
}

template <class BlockedIn, class BlockedOut>
void device_copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out) {
    if (pluto::devices() == 0) {
        return host_copy_blocked_to_blocked_mdspan(blocked_in, blocked_out);
    }
    ATLAS_NOTIMPLEMENTED;
}
#else
// Implementation in transpositions_on_device.hic
template <class Nonblocked, class Blocked>
void device_copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked);

template <class Blocked, class Nonblocked>
void device_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked);

template <class BlockedIn, class BlockedOut>
void device_copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out);
#endif

//

template <class Nonblocked, class Blocked>
void copy_nonblocked_to_blocked_mdspan(const Nonblocked nonblocked, Blocked blocked, bool on_device) {
    ATLAS_TRACE("copy_nonblocked_to_blocked_mdspan "+std::string(on_device?"[device]":"[host]"));
    if (on_device) {
        device_copy_nonblocked_to_blocked_mdspan(nonblocked, blocked);
    }
    else {
        host_copy_nonblocked_to_blocked_mdspan(nonblocked, blocked);
    }
}

template <class Blocked, class Nonblocked>
void copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked, bool on_device) {
    ATLAS_TRACE("copy_blocked_to_nonblocked_mdspan "+std::string(on_device?"[device]":"[host]"));
    if (on_device) {
        device_copy_blocked_to_nonblocked_mdspan(blocked, nonblocked);
    }
    else {
        host_copy_blocked_to_nonblocked_mdspan(blocked, nonblocked);
    }
}

template <class BlockedIn, class BlockedOut>
void copy_blocked_to_blocked_mdspan(const BlockedIn blocked_in, BlockedOut blocked_out, bool on_device) {
    ATLAS_TRACE("copy_blocked_to_blocked_mdspan "+std::string(on_device?"[device]":"[host]"));
    if (on_device) {
        device_copy_blocked_to_blocked_mdspan(blocked_in, blocked_out);
    }
    else {
        host_copy_blocked_to_blocked_mdspan(blocked_in, blocked_out);
    }
}

// Routines that delegate

template <class ValueType>
void copy_blocked_to_nonblocked_T(const array::Array& blocked, array::Array& nonblocked, bool on_device) {
    ATLAS_ASSERT(nonblocked.rank() == blocked.rank()-1);
    if (blocked.rank()==4) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType, 4>(blocked)    : array::make_host_view<ValueType, 4>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType, 3>(nonblocked) : array::make_host_view<ValueType, 3>(nonblocked);
        if(USE_MDSPAN) {
            copy_blocked_to_nonblocked_mdspan(blocked_v.as_mdspan(), nonblocked_v.as_mdspan(), on_device);
        }
        else {
            copy_blocked_to_nonblocked_mdspan(blocked_v, nonblocked_v, on_device);
        }
    }
    else if (blocked.rank()==3) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType, 3>(blocked)    : array::make_host_view<ValueType, 3>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType, 2>(nonblocked) : array::make_host_view<ValueType, 2>(nonblocked);
        if(USE_MDSPAN) {
            copy_blocked_to_nonblocked_mdspan(blocked_v.as_mdspan(), nonblocked_v.as_mdspan(), on_device);
        }
        else {
            copy_blocked_to_nonblocked_mdspan(blocked_v, nonblocked_v, on_device);
        }
    }
    else if (blocked.rank()==2) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType, 2>(blocked)    : array::make_host_view<ValueType, 2>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType, 1>(nonblocked) : array::make_host_view<ValueType, 1>(nonblocked);
        if(USE_MDSPAN) {
            copy_blocked_to_nonblocked_mdspan(blocked_v.as_mdspan(), nonblocked_v.as_mdspan(), on_device);
        }
        else {
            copy_blocked_to_nonblocked_mdspan(blocked_v, nonblocked_v, on_device);
        }
    }
    else {
        ATLAS_THROW_EXCEPTION("transposition not implemented");
    }
}

template <class ValueType>
void copy_nonblocked_to_blocked_T(const array::Array& nonblocked, array::Array& blocked, bool on_device) {
    ATLAS_ASSERT(nonblocked.rank() == blocked.rank()-1);
    if (blocked.rank()==4) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType,4>(blocked)    : array::make_host_view<ValueType, 4>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType,3>(nonblocked) : array::make_host_view<ValueType, 3>(nonblocked);
        if(USE_MDSPAN) {
            copy_nonblocked_to_blocked_mdspan(nonblocked_v.as_mdspan(), blocked_v.as_mdspan(), on_device);
        }
        else {
            copy_nonblocked_to_blocked_mdspan(nonblocked_v, blocked_v, on_device);
        }
    }
    else if (blocked.rank()==3) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType,3>(blocked)    : array::make_host_view<ValueType, 3>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType,2>(nonblocked) : array::make_host_view<ValueType, 2>(nonblocked);
        if(USE_MDSPAN) {
            copy_nonblocked_to_blocked_mdspan(nonblocked_v.as_mdspan(), blocked_v.as_mdspan(), on_device);
        }
        else {
            copy_nonblocked_to_blocked_mdspan(nonblocked_v, blocked_v, on_device);
        }
    }
    else if (blocked.rank()==2) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType,2>(blocked)    : array::make_host_view<ValueType, 2>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType,1>(nonblocked) : array::make_host_view<ValueType, 1>(nonblocked);
        if(USE_MDSPAN) {
            copy_nonblocked_to_blocked_mdspan(nonblocked_v.as_mdspan(), blocked_v.as_mdspan(), on_device);
        }
        else {
            copy_nonblocked_to_blocked_mdspan(nonblocked_v, blocked_v, on_device);
        }
    }
    else {
        ATLAS_THROW_EXCEPTION("transposition not implemented");
    }
}

template <class ValueType>
void copy_blocked_to_blocked_T(const array::Array& blocked_in, array::Array& blocked_out, bool on_device) {
    ATLAS_ASSERT(blocked_in.rank() == blocked_out.rank());
    if (blocked_in.rank()==4) {
        auto blocked_in_v  = on_device ? array::make_device_view<ValueType, 4>(blocked_in)  : array::make_host_view<ValueType, 4>(blocked_in);
        auto blocked_out_v = on_device ? array::make_device_view<ValueType, 4>(blocked_out) : array::make_host_view<ValueType, 4>(blocked_out);
        if (USE_MDSPAN) {
            copy_blocked_to_blocked_mdspan(blocked_in_v.as_mdspan(), blocked_out_v.as_mdspan(), on_device);
        }
        else {
            copy_blocked_to_blocked_mdspan(blocked_in_v, blocked_out_v, on_device);
        }
    }
    else if (blocked_in.rank()==3) {
        auto blocked_in_v  = on_device ? array::make_device_view<ValueType, 3>(blocked_in)  : array::make_host_view<ValueType, 3>(blocked_in);
        auto blocked_out_v = on_device ? array::make_device_view<ValueType, 3>(blocked_out) : array::make_host_view<ValueType, 3>(blocked_out);
        if (USE_MDSPAN) {
            copy_blocked_to_blocked_mdspan(blocked_in_v.as_mdspan(), blocked_out_v.as_mdspan(), on_device);
        }
        else {
            copy_blocked_to_blocked_mdspan(blocked_in_v, blocked_out_v, on_device);
        }
    }
    else if (blocked_in.rank()==2) {
        auto blocked_in_v  = on_device ? array::make_device_view<ValueType, 2>(blocked_in)  : array::make_host_view<ValueType, 2>(blocked_in);
        auto blocked_out_v = on_device ? array::make_device_view<ValueType, 2>(blocked_out) : array::make_host_view<ValueType, 2>(blocked_out);
        if (USE_MDSPAN) {
            copy_blocked_to_blocked_mdspan(blocked_in_v.as_mdspan(), blocked_out_v.as_mdspan(), on_device);
        }
        else {
            copy_blocked_to_blocked_mdspan(blocked_in_v, blocked_out_v, on_device);
        }
    }
    else {
        ATLAS_THROW_EXCEPTION("transposition not implemented");
    }
}


void copy_nonblocked_to_blocked(const array::Array& nonblocked, array::Array& blocked, bool on_device) {
    ATLAS_ASSERT(blocked.datatype() == nonblocked.datatype());
    switch (blocked.datatype().kind()) {
        case array::DataType::kind<int>()    : return copy_nonblocked_to_blocked_T<int>(nonblocked, blocked, on_device);
        case array::DataType::kind<long>()   : return copy_nonblocked_to_blocked_T<long>(nonblocked, blocked, on_device);
        case array::DataType::kind<float>()  : return copy_nonblocked_to_blocked_T<float>(nonblocked, blocked, on_device);
        case array::DataType::kind<double>() : return copy_nonblocked_to_blocked_T<double>(nonblocked, blocked, on_device);
        default: throw_Exception("datatype not supported", Here());
    }
}

void copy_blocked_to_nonblocked(const array::Array& blocked, array::Array& nonblocked, bool on_device) {
    ATLAS_ASSERT(blocked.datatype() == nonblocked.datatype());
    switch (blocked.datatype().kind()) {
        case array::DataType::kind<int>()    : return copy_blocked_to_nonblocked_T<int>(blocked, nonblocked, on_device);
        case array::DataType::kind<long>()   : return copy_blocked_to_nonblocked_T<long>(blocked, nonblocked, on_device);
        case array::DataType::kind<float>()  : return copy_blocked_to_nonblocked_T<float>(blocked, nonblocked, on_device);
        case array::DataType::kind<double>() : return copy_blocked_to_nonblocked_T<double>(blocked, nonblocked, on_device);
        default: throw_Exception("datatype not supported", Here());
    }
}

void copy_blocked_to_blocked(const array::Array& blocked_in, array::Array& blocked_out, bool on_device) {
    ATLAS_ASSERT(blocked_in.datatype() == blocked_out.datatype());
    switch (blocked_in.datatype().kind()) {
        case array::DataType::kind<int>()    : return copy_blocked_to_blocked_T<int>(blocked_in, blocked_out, on_device);
        case array::DataType::kind<long>()   : return copy_blocked_to_blocked_T<long>(blocked_in, blocked_out, on_device);
        case array::DataType::kind<float>()  : return copy_blocked_to_blocked_T<float>(blocked_in, blocked_out, on_device);
        case array::DataType::kind<double>() : return copy_blocked_to_blocked_T<double>(blocked_in, blocked_out, on_device);
        default: throw_Exception("datatype not supported", Here());
    }
}

void copy_blocked_to_blocked(const FieldSet& blocked_fields_in, FieldSet& blocked_fields_out, bool on_device) {
    ATLAS_ASSERT(blocked_fields_in.size() == blocked_fields_out.size());
    for( int i=0; i<blocked_fields_in.size(); ++i) {
        copy_blocked_to_blocked(blocked_fields_in[i], blocked_fields_out[i], on_device);
    }
}

void copy_blocked_to_nonblocked(const FieldSet& blocked_fields_in, FieldSet& nonblocked_fields_out, bool on_device) {
    ATLAS_ASSERT(blocked_fields_in.size() == nonblocked_fields_out.size());
    for( int i=0; i<blocked_fields_in.size(); ++i) {
        copy_blocked_to_nonblocked(blocked_fields_in[i], nonblocked_fields_out[i], on_device);
    }
}

void copy_nonblocked_to_blocked(const FieldSet& nonblocked_fields_in, FieldSet& blocked_fields_out, bool on_device) {
    ATLAS_ASSERT(nonblocked_fields_in.size() == blocked_fields_out.size());
    for( int i=0; i<blocked_fields_out.size(); ++i) {
        copy_nonblocked_to_blocked(nonblocked_fields_in[i], blocked_fields_out[i], on_device);
    }
}

}

