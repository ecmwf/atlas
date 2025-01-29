
/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "pointer_info.h"

#include <tuple> // std::ignore

#include "hic/hic.h"

#include "pluto/pluto_config.h"
#include "pluto/runtime.h"


#define LOG PLUTO_DEBUGGING

namespace pluto {

#ifdef HIC_CALL_OR_RETURN
#undef HIC_CALL_OR_RETURN
#endif

#define HIC_CALL_OR_RETURN( code, value ) \
    if( code != hicSuccess ) { \
        std::ignore = hicGetLastError(); \
        return value; \
    }

bool is_device_accessible(const void* ptr) {
    constexpr bool default_value = false;
    if constexpr(PLUTO_HAVE_HIC) {
        if (not devices()) {
            return default_value;
        }
        hicPointerAttributes attr;
        HIC_CALL_OR_RETURN( hicPointerGetAttributes(&attr, ptr), default_value );
        return ptr == attr.devicePointer;
    }
    return default_value;
}
bool is_host_accessible(const void* ptr) {
    constexpr bool default_value = true;
    if constexpr(PLUTO_HAVE_HIC) {
        if (not devices()) {
            return default_value;
        }
        hicPointerAttributes attr;
        HIC_CALL_OR_RETURN( hicPointerGetAttributes(&attr, ptr), default_value );
        return ptr == attr.hostPointer;
    }
    return default_value;
}

bool is_managed(const void* ptr) {
    constexpr bool default_value = false;
    if constexpr(PLUTO_HAVE_HIC) {
        if (not devices()) {
            return default_value;
        }
        hicPointerAttributes attr;
        HIC_CALL_OR_RETURN( hicPointerGetAttributes(&attr, ptr), default_value );
        return attr.type == hicMemoryTypeManaged;
    }
    return default_value;
}

bool is_pinned(const void* ptr) {
    constexpr bool default_value = false;
    if constexpr(PLUTO_HAVE_HIC) {
        if (not devices()) {
            return default_value;
        }
        hicPointerAttributes attr;
        HIC_CALL_OR_RETURN( hicPointerGetAttributes(&attr, ptr), default_value );
        return (attr.type != hicMemoryTypeManaged
            &&  attr.type != hicMemoryTypeUnregistered
            &&  attr.hostPointer != nullptr
            &&  attr.devicePointer != nullptr);
    }
    return default_value;
}

bool is_host(const void* ptr) {
    constexpr bool default_value = true;
    if constexpr(PLUTO_HAVE_HIC) {
        if (not devices()) {
            return default_value;
        }
        hicPointerAttributes attr;
        HIC_CALL_OR_RETURN( hicPointerGetAttributes(&attr, ptr), default_value );
        return attr.type == hicMemoryTypeUnregistered || attr.type == hicMemoryTypeHost;
    }
    return default_value;
}


bool is_device(const void* ptr) {
    constexpr bool default_value = false;
    if constexpr(PLUTO_HAVE_HIC) {
        if (not devices()) {
            return default_value;
        }
        hicPointerAttributes attr;
        HIC_CALL_OR_RETURN( hicPointerGetAttributes(&attr, ptr), default_value);
        return attr.type == hicMemoryTypeDevice;
    }
    return default_value;
}


void* get_registered_device_pointer(const void* host_ptr) {
    if constexpr(PLUTO_HAVE_HIC) {
        if (not devices()) {
            return const_cast<void*>(host_ptr);
        }
        void* device_ptr;
        auto code = hicHostGetDevicePointer(&device_ptr, const_cast<void*>(host_ptr), 0);
        std::ignore = hicGetLastError(); // sets stored error code back to success
        if (code != hicSuccess) {
            return const_cast<void*>(host_ptr);
        }
        return device_ptr;
    }
    return const_cast<void*>(host_ptr);
}

}
