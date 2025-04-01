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

#include "pluto/pluto_config.h"

#include "pluto/memory.h"
#include "pluto/memory_resource.h"
#include "pluto/memory_resource/AsyncMemoryResourceAdaptor.h"
#include "pluto/memory_resource/DeviceMemoryResource.h"
#include "pluto/memory_resource/HostMemoryResource.h"
#include "pluto/memory_resource/ManagedMemoryResource.h"
#include "pluto/memory_resource/MemoryPoolResource.h"
#include "pluto/memory_resource/MemoryResourceAdaptor.h"
#include "pluto/memory_resource/PinnedMemoryResource.h"
#include "pluto/memory_resource/TraceMemoryResource.h"

#include "pluto/device/MemoryResource.h"
#include "pluto/device/allocator.h"
#include "pluto/host/MemoryResource.h"
#include "pluto/host/allocator.h"

#include "pluto/copy.h"
#include "pluto/event.h"
#include "pluto/memcpy.h"
#include "pluto/prefetch.h"
#include "pluto/stream.h"
#include "pluto/wait.h"

#include "pluto/alignment.h"
#include "pluto/pointer_info.h"
#include "pluto/runtime.h"
#include "pluto/scope.h"

namespace pluto {
inline void release() {
    if (trace::enabled()) {
        trace::out << "pluto::release()" << std::endl;
    }
    pluto::host_pool_resource()->release();
    pluto::pinned_pool_resource()->release();
    pluto::device_pool_resource()->release();
    pluto::managed_pool_resource()->release();
}
}
