/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "MemoryResource.h"

#include <cstdlib>

#include "pluto/memory_resource.h"
#include "pluto/memory_resource/DeviceMemoryResource.h"
#include "pluto/pluto_config.h"

#define LOG PLUTO_DEBUGGING

// ---------------------------------------------------------------------------------------------------------

namespace pluto::device {

// ---------------------------------------------------------------------------------------------------------

memory_resource* default_{nullptr};

memory_resource* get_default_resource() {
    if (default_ == nullptr) {
        default_ = device_resource();
        if (const auto* env = std::getenv("PLUTO_DEVICE_MEMORY_RESOURCE")) {
            default_ = get_registered_resource(env);
        }
    }
    return default_;
}

memory_resource* set_default_resource(memory_resource* mr) {
    auto* previous = default_;
    default_ = mr;
    return previous;
}

memory_resource* set_default_resource(std::string_view name) {
    return pluto::device::set_default_resource(get_registered_resource(name));
}

// ---------------------------------------------------------------------------------------------------------

scoped_default_resource::scoped_default_resource(std::string_view name):
    scoped_default_resource(get_registered_resource(name)) {}

// ---------------------------------------------------------------------------------------------------------


}  // namespace pluto::device
