/*
* (C) Copyright 2024 ECMWF.
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation nor
* does it submit to any jurisdiction.
*/

#include "RegisterPointerInfo.h"

#include <map>
#include <string>

namespace atlas {
namespace util {

//------------------------------------------------------------------------------

static std::map<const void*, std::string, std::less<>> map_pointer_name;

void register_pointer_name(const void* pointer, std::string_view name) {
    map_pointer_name.insert_or_assign(pointer, name);
}

void unregister_pointer_name(const void* pointer) {
    map_pointer_name.erase(pointer);
}

std::string_view registered_pointer_name(const void* pointer) {
    auto it = map_pointer_name.find(pointer);
    if( it != map_pointer_name.end() ) {
        return it->second;
    }
    return {};
}

//------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
