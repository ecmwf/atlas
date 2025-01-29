
/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "prefetch.h"

#include "hic/hic.h"
#include "pluto/pluto_config.h"

#include "pluto/stream.h"

namespace pluto {

void prefetch_host_to_device(const void* managed_ptr, std::size_t bytes) {
    if constexpr(PLUTO_HAVE_HIC) {
        HIC_CALL( hicMemPrefetchAsync(managed_ptr, bytes, 0 /*device id*/) );
    }
}

void prefetch_host_to_device(const void* managed_ptr, std::size_t bytes, const stream& s) {
    if constexpr(PLUTO_HAVE_HIC) {
        HIC_CALL( hicMemPrefetchAsync(managed_ptr, bytes, 0 /*device id*/, s.value<hicStream_t>() ) );
    }
}

void prefetch_device_to_host(const void* managed_ptr, std::size_t bytes) {
    if constexpr(PLUTO_HAVE_HIC) {
        HIC_CALL( hicMemPrefetchAsync(managed_ptr, bytes, hicCpuDeviceId) );
    }
}

void prefetch_device_to_host(const void* managed_ptr, std::size_t bytes, const stream& s) {
    if constexpr(PLUTO_HAVE_HIC) {
        HIC_CALL( hicMemPrefetchAsync(managed_ptr, bytes, hicCpuDeviceId, s.value<hicStream_t>()) );
    }
}

}
