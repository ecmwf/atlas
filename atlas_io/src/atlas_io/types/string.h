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

#include <string>

#include "atlas_io/Data.h"
#include "atlas_io/Metadata.h"
#include "atlas_io/atlas_compat.h"
#include "atlas_io/detail/Assert.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

inline size_t encode_metadata(const std::string& value, atlas::io::Metadata& out) {
    out.set("type", "string");
    out.set("value", value);
    return 0;
}

inline void encode_data(const std::string&, atlas::io::Data&) {}

inline void decode(const atlas::io::Metadata& metadata, const atlas::io::Data&, std::string& value) {
    ATLAS_IO_ASSERT(metadata.getString("type") == "string");
    metadata.get("value", value);
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
