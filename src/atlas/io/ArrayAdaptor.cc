/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "ArrayAdaptor.h"

#include <cstring>  // memcpy
#include <sstream>

#include "atlas/array/Array.h"

namespace atlas {
namespace array {

//---------------------------------------------------------------------------------------------------------------------

void interprete(const atlas::array::Array& a, atlas::io::ArrayReference& out) {
    out = io::ArrayReference(a.data(), atlas::io::DataType(a.datatype().str()), a.shape());
}

//---------------------------------------------------------------------------------------------------------------------

void decode(const atlas::io::Metadata& metadata, const atlas::io::Data& data, atlas::array::Array& out) {
    atlas::io::ArrayMetadata array(metadata);

    if (array.datatype().str() != out.datatype().str()) {
        std::stringstream err;
        err << "Could not decode " << metadata.json() << " into Array with datatype " << out.datatype().str() << "."
            << "Incompatible datatype!";
        throw atlas::io::Exception(err.str(), Here());
    }
    if (array.rank() != out.rank()) {
        std::stringstream err;
        err << "Could not decode " << metadata.json() << " into Array with rank " << out.rank() << "."
            << "Incompatible rank!";
        throw atlas::io::Exception(err.str(), Here());
    }

    out.resize(array.shape());

    ATLAS_IO_ASSERT(out.contiguous());

    ::memcpy(out.data(), data, data.size());
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
