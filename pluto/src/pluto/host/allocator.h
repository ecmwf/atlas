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

#include "pluto/host/MemoryResource.h"

namespace pluto::host {

// --------------------------------------------------------------------------------------------------------

template <typename T>
class allocator : public pluto::allocator<T> {
public:
    using pluto::allocator<T>::allocator;
    allocator(): pluto::allocator<T>::allocator(host::get_default_resource()) {}
};

// --------------------------------------------------------------------------------------------------------

}  // namespace pluto::host
