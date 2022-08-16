/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Metadata.h"

#include <limits>
#include <ostream>
#include <sstream>

#include "eckit/log/JSON.h"

#include "atlas_io/atlas_compat.h"
#include "atlas_io/detail/Assert.h"

#include "atlas_io/Exceptions.h"
#include "atlas_io/types/array/ArrayReference.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

size_t uncompressed_size(const atlas::io::Metadata& m) {
    if (m.has("data.size")) {
        return m.getUnsigned("data.size");
    }
    else if (m.has("type")) {
        if (m.getString("type") == "array") {
            atlas::io::ArrayMetadata array(m);
            return array.bytes();
        }
    }
    std::stringstream err;
    err << "Could not compute uncompressed data size from metadata \n";
    write(m, err);
    throw Exception(err.str());
}

//---------------------------------------------------------------------------------------------------------------------

void write(const atlas::io::Metadata& metadata, std::ostream& out) {
    eckit::JSON js(out, eckit::JSON::Formatting::indent(4));
    js << metadata;
}

void write(const atlas::io::Metadata& metadata, atlas::io::Stream& out) {
    std::stringstream ss;
    write(metadata, ss);
    std::string s = ss.str();
    out.write(s.data(), s.size());
}

//---------------------------------------------------------------------------------------------------------------------

void Metadata::link(Metadata&& linked) {
    std::string initial_link = link();
    ATLAS_IO_ASSERT(initial_link.size());

    data   = std::move(linked.data);
    record = std::move(linked.record);

    set(linked);

    // Set link to initial_link, in case the link is itself another link
    set("link", initial_link);
}

std::string Metadata::json() const {
    std::stringstream s;
    eckit::JSON js(s, eckit::JSON::Formatting::compact());
    js << *this;
    return s.str();
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
