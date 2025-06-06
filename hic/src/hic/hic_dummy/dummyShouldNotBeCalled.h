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

#include <stdexcept>
#include <string>

namespace {

[[noreturn]] void dummyShouldNotBeCalled(const char* symbol) {
    throw std::runtime_error(std::string(symbol) + " is using the dummy backend and should not be called");
}

}  // namespace