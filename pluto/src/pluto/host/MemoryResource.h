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

#include <string_view>
#include <cassert>

#include "pluto/memory_resource/memory_resource.h"

namespace pluto::host {

// --------------------------------------------------------------------------------------------------------

inline void set_default_resource(memory_resource& mr) {
    pluto::set_default_resource(&mr);
}

inline void set_default_resource(memory_resource* mr) {
    pluto::set_default_resource(mr);
}

inline void set_default_resource(std::string_view name) {
    pluto::set_default_resource(get_registered_resource(name));
}

inline memory_resource* get_default_resource() {
    pluto::init();
    return pluto::get_default_resource();
}

class [[nodiscard]] DefaultResource {
public:
    DefaultResource(std::string_view name) : DefaultResource(get_registered_resource(name)) {}

    DefaultResource(memory_resource* mr) :
        saved_(pluto::host::get_default_resource()) {
        pluto::host::set_default_resource(mr);
    }
    ~DefaultResource() {
        pluto::host::set_default_resource(saved_);
    }
private:
    memory_resource* saved_;
};

// --------------------------------------------------------------------------------------------------------

}

