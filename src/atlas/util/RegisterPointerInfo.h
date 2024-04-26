/*
* (C) Copyright 2013 ECMWF.
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation nor
* does it submit to any jurisdiction.
*/

#pragma once

#include <string_view>

namespace atlas {
namespace util {

//------------------------------------------------------------------------------

void register_pointer_name(const void* pointer, std::string_view name);
void unregister_pointer_name(const void* pointer);
std::string_view registered_pointer_name(const void* pointer);

//------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
