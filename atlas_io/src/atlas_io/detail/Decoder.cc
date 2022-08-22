/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Decoder.h"

#include "atlas_io/Trace.h"

namespace atlas {
namespace io {

void decode(const atlas::io::Metadata& metadata, const atlas::io::Data& data, Decoder& decoder) {
    ATLAS_IO_TRACE("decode");
    decoder.self_->decode_(metadata, data);
}

void decode(const atlas::io::Metadata& metadata, const atlas::io::Data& data, Decoder&& decoder) {
    ATLAS_IO_TRACE_SCOPE("decode");
    decoder.self_->decode_(metadata, data);
}


}  // namespace io
}  // namespace atlas
