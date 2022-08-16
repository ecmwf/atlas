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

#include "atlas_io/Record.h"
#include "atlas_io/Stream.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

class Session {
public:
    Session();

    ~Session();

    static bool active();

    static Record record(const std::string& path, size_t offset);

    static Record record(Stream, size_t offset);

    static void store(Stream stream);
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
