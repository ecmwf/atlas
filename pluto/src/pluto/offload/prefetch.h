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

#include <cstddef>

#include "pluto/offload/Stream.h"

namespace pluto {

void prefetch_managed_device(const void* managed_ptr, std::size_t bytes);
void prefetch_managed_device(const void* managed_ptr, std::size_t bytes, const Stream&);

void prefetch_managed_host(const void* managed_ptr, std::size_t bytes);
void prefetch_managed_host(const void* managed_ptr, std::size_t bytes, const Stream&);

}
