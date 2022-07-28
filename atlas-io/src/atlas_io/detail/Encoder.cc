/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Encoder.h"

#include "atlas_io/Trace.h"

namespace atlas {
namespace io {

size_t encode_metadata(const Encoder& encoder, atlas::io::Metadata& metadata) {
    ATLAS_IO_TRACE();
    ASSERT(encoder);
    return encoder.self_->encode_metadata_(metadata);
}

void encode_data(const Encoder& encoder, atlas::io::Data& out) {
    ATLAS_IO_TRACE();
    ASSERT(encoder);
    encoder.self_->encode_data_(out);
}

}  // namespace io
}  // namespace atlas
