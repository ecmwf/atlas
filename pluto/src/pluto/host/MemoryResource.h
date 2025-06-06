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

#include <cassert>
#include <string_view>

#include "pluto/memory_resource.h"

namespace pluto::host {

// --------------------------------------------------------------------------------------------------------

memory_resource* get_default_resource();
memory_resource* set_default_resource(memory_resource*);
memory_resource* set_default_resource(std::string_view name);


class [[nodiscard]] scoped_default_resource {
public:
    scoped_default_resource(std::string_view name);

    scoped_default_resource(memory_resource* mr): saved_(get_default_resource()) { host::set_default_resource(mr); }

    ~scoped_default_resource() { host::set_default_resource(saved_); }

private:
    memory_resource* saved_;
};

// --------------------------------------------------------------------------------------------------------

}  // namespace pluto::host
