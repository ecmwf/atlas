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
#include <cstdlib>

namespace atlas::acc {

enum class CompilerId {
    unknown,
    nvidia,
    cray,
};

int devices();
void map(void* host_data, void* device_data, std::size_t bytes);
void unmap(void* host_data);
bool is_present(void* host_data, std::size_t bytes);
void* deviceptr(void* host_data);
CompilerId compiler_id();

}

