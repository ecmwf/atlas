/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstdint>
#include <string>

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

class Base64 {
public:
    static std::string encode(const void* data, size_t len);
    static std::string decode(const void* data, size_t len);

    template <typename T>
    static std::string encode(const T& value) {
        return encode(&value, sizeof(value));
    }

    template <typename T>
    static T decode(const std::string& in) {
        std::string decoded = decode(in.data(), in.size());
        return *reinterpret_cast<const T*>(decoded.data());
    }
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
