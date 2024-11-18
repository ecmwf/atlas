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

#include "pluto/memory_resource/memory_resource.h"
#include "pluto/memory_resource/MemoryResourceAdaptor.h"
#include "pluto/memory_resource/TraceMemoryResource.h"
#include "pluto/memory_resource/MemoryPoolResource.h"
#include "pluto/memory_resource/PinnedMemoryResource.h"
#include "pluto/memory_resource/ManagedMemoryResource.h"
#include "pluto/memory_resource/DeviceMemoryResource.h"
#include "pluto/memory_resource/StreamMemoryResourceAdaptor.h"

#include "pluto/host/allocator.h"
#include "pluto/host/MemoryResource.h"
#include "pluto/host/unique_ptr.h"
#include "pluto/host/make_copy.h"
#include "pluto/host/vector.h"
#include "pluto/device/MemoryResource.h"
#include "pluto/device/allocator.h"
#include "pluto/device/unique_ptr.h"
#include "pluto/device/make_copy.h"
#include "pluto/device/vector.h"

#include "pluto/offload/copy.h"
#include "pluto/offload/Event.h"
#include "pluto/offload/memcpy.h"
#include "pluto/offload/prefetch.h"
#include "pluto/offload/Stream.h"
#include "pluto/offload/wait.h"

#include "pluto/util/Alignment.h"
#include "pluto/util/PointerInfo.h"
#include "pluto/util/Runtime.h"
#include "pluto/util/Scope.h"
